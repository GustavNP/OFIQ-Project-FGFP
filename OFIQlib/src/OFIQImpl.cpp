/**
 * @file OFIQImpl.cpp
 *
 * @copyright Copyright (c) 2024  Federal Office for Information Security, Germany
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * @author OFIQ development team
 */




/* ========== FINE-GRAINED FACE PARSING ============
    This file has been modified from the original OFIQ project.
    It now includes a Fine-Grained Face Parsing algorithm, dividing the face parsing mask into smaller regions than the existing algorithm does.
    
    Author: Gustav Nilsson Pedersen - s174562@student.dtu.dk
*/


#include "Configuration.h"
#include "Executor.h"
#include "ofiq_lib_impl.h"
#include "OFIQError.h"
#include "FaceMeasures.h"
#include "utils.h"
#include "image_io.h"
#include <chrono>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>
#include <map>
using hrclock = std::chrono::high_resolution_clock;

using namespace std;
using namespace OFIQ;
using namespace OFIQ_LIB;
using namespace OFIQ_LIB::modules::measures;


OFIQImpl::OFIQImpl():m_emptySession({this->dummyImage, this->dummyAssement}) {}

ReturnStatus OFIQImpl::initialize(const std::string& configDir, const std::string& configFilename)
{
    try
    {
        this->config = std::make_unique<Configuration>(configDir, configFilename);
        CreateNetworks();
        m_executorPtr = CreateExecutor(m_emptySession);
    }
    catch (const OFIQError & ex)
    {
        return {ex.whatCode(), ex.what()};
    }
    catch (const std::exception & ex)
    {
        return {ReturnCode::UnknownError, ex.what()};
    }

    return ReturnStatus(ReturnCode::Success);
}

ReturnStatus OFIQImpl::scalarQuality(const OFIQ::Image& face, double& quality, const std::string& inputFile)
{
    FaceImageQualityAssessment assessments;

    if (auto result = vectorQuality(face, assessments, inputFile);
        result.code != ReturnCode::Success)
        return result;

    if(assessments.qAssessments.find(QualityMeasure::UnifiedQualityScore) != assessments.qAssessments.end())
        quality = assessments.qAssessments[QualityMeasure::UnifiedQualityScore].scalar;
    else
    {
        // just as an example - the scalarQuality is an average of all valid scalar measurements
        double sumScalars = 0;
        int numScalars = 0;
        for (auto const& [measureName, aResult] : assessments.qAssessments)
        {
            if (aResult.scalar != -1)
            {
                sumScalars += aResult.scalar;
                numScalars++;
            }
        }
        quality = numScalars != 0 ? sumScalars / numScalars : 0;
    }

    return ReturnStatus(ReturnCode::Success);
}

void OFIQImpl::performPreprocessing(Session& session, const std::string& inputFile)
{
    std::chrono::time_point<hrclock> tic;

    log("\t1. detectFaces ");
    tic = hrclock::now();

    std::vector<OFIQ::BoundingBox> faces = networks->faceDetector->detectFaces(session);
    if (faces.empty())
    {
        log("\n\tNo faces were detected, abort preprocessing\n");
        throw OFIQError(ReturnCode::FaceDetectionError, "No faces were detected");
    }
    log(std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            hrclock::now() - tic).count()) + std::string(" ms "));

    session.setDetectedFaces(faces);
    log("2. estimatePose ");
    tic = hrclock::now();

    session.setPose(networks->poseEstimator->estimatePose(session));

    log(std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            hrclock::now() - tic).count()) + std::string(" ms "));

    log("3. extractLandmarks ");
    tic = hrclock::now();

    session.setLandmarks(networks->landmarkExtractor->extractLandmarks(session));

    log(std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            hrclock::now() - tic).count()) + std::string(" ms "));

    log("4. alignFaceImage ");
    tic = hrclock::now();
    // aligned face requires the landmarks of the face thus it must come after the landmark extraction.
    alignFaceImage(session);
    log(std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            hrclock::now() - tic).count()) + std::string(" ms "));

    log("5. getSegmentationMask ");
    tic = hrclock::now();
    // segmentation results for face_parsing
    session.setFaceParsingImage(OFIQ_LIB::copyToCvImage(
        networks->segmentationExtractor->GetMask(
            session,
            OFIQ_LIB::modules::segmentations::SegmentClassLabels::face),
        true));
    log(std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            hrclock::now() - tic).count()) + std::string(" ms "));
    
    //  ==========  FINE-GRAINED FACE PARSING ALGORITHM ========
    fineGrainedFaceParsing(session, inputFile);

    log("6. getFaceOcclusionMask ");
    tic = hrclock::now();
    session.setFaceOcclusionSegmentationImage(OFIQ_LIB::copyToCvImage(
        networks->faceOcclusionExtractor->GetMask(
            session,
            OFIQ_LIB::modules::segmentations::SegmentClassLabels::face),
        true));
    log(std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            hrclock::now() - tic).count()) + std::string(" ms "));

    static const std::string alphaParamPath = "params.measures.FaceRegion.alpha";
    double alpha = 0.0f;
    try
    {
        alpha = this->config->GetNumber(alphaParamPath);
    }
    catch(...)
    {
        alpha = 0.0f;
    }

    log("7. getAlignedFaceMask ");
    tic = hrclock::now();

    session.setAlignedFaceLandmarkedRegion(
         OFIQ_LIB::modules::landmarks::FaceMeasures::GetFaceMask(
            session.getAlignedFaceLandmarks(),
            session.getAlignedFace().rows,
            session.getAlignedFace().cols,
            (float)alpha
         )
    );
    log(std::to_string(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            hrclock::now() - tic).count()) + std::string(" ms "));

    log("\npreprocessing finished\n");
}

void OFIQImpl::alignFaceImage(Session& session) {
    auto landmarks = session.getLandmarks();
    OFIQ::FaceLandmarks alignedFaceLandmarks;
    alignedFaceLandmarks.type = landmarks.type;
    cv::Mat transformationMatrix;
    cv::Mat alignedBGRimage = alignImage(session.image(), landmarks, alignedFaceLandmarks, transformationMatrix);

    session.setAlignedFace(alignedBGRimage);
    session.setAlignedFaceLandmarks(alignedFaceLandmarks);
    session.setAlignedFaceTransformationMatrix(transformationMatrix);

}


// === Part of the Fine-Grained Face Parsing algorithm ===
cv::Mat OFIQImpl::getPolygonMask(OFIQ::Landmarks polygonPoints, int maskSizeOneDirection, int pixelValue) {

    int noOfPoints = static_cast<int>(polygonPoints.size());

    std::vector<cv::Point> landmarkPolygonPoints;

    for (int i = 0; i < noOfPoints; i++)
    {
        cv::Point point = cv::Point(polygonPoints[i].x, polygonPoints[i].y);
        landmarkPolygonPoints.push_back(point);
    }

    std::vector<vector<cv::Point>> vectorOfVectorWithPoints; // cv::fillPoly accepts vector of vectors with points for polygons (so multiple polygons)
    vectorOfVectorWithPoints.push_back(landmarkPolygonPoints);

    cv::Mat mask = cv::Mat::zeros(maskSizeOneDirection, maskSizeOneDirection, CV_8U);

    cv::fillPoly(mask, vectorOfVectorWithPoints, cv::Scalar(pixelValue));


    return mask;
}

// ======== Fine-Grained Face Parsing algorithm ========
void OFIQImpl::fineGrainedFaceParsing(Session& session, const string& inputFile) {
    auto landmarks = session.getAlignedFaceLandmarks();

    cv::Mat faceParsingImage = session.getFaceParsingImage().clone();
    cv::resize(faceParsingImage, faceParsingImage, cv::Size(200, 200), 0, 0, cv::INTER_NEAREST);

    float scalingFactor = 200.0f / (session.getAlignedFace().rows - 60.0f);

    OFIQ::Landmarks scaledLandmarks;
    // for each landmark subtract 30 from the x-value
    // multiply the x- and y-value in each landmark with the scaling factor
    for (int i = 0; i < landmarks.landmarks.size(); i++)
    {
        auto newX = static_cast<int>(std::round((static_cast<float>(landmarks.landmarks[i].x) - 30) * scalingFactor));
        auto newY = static_cast<int>(std::round(static_cast<float>(landmarks.landmarks[i].y) * scalingFactor));

        scaledLandmarks.emplace_back(
            OFIQ::LandmarkPoint(static_cast<uint16_t>(newX), static_cast<uint16_t>(newY)));
    }


    std::vector<int> rightOrbitalRegionLandmarkIds{ 33, 41, 40, 39, 38, 51, 64, 63, 62, 61, 60 };
    std::vector<int> leftOrbitalRegionLandmarkIds{ 46, 47, 48, 49, 50, 51, 68, 69, 70, 71, 72 };
    std::vector<int> nasalRegionLandmarkIds{ 51, 64, 76, 77, 78, 79, 80, 81, 82, 68 };
    std::vector<int> mentalRegionLandmarkIds{ 82, 83, 84, 85, 86, 87, 76, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };
    std::vector<int> rightBuccalRegionLandmarkIds{ 5, 6, 7, 8, 9, 10, 11, 76, 55 };
    std::vector<int> leftBuccalRegionLandmarkIds{ 27, 26, 25, 24, 23, 22, 21, 82, 59 };
    std::vector<int> rightZygoInfraParoRegionLandmarkIds{ 0, 1, 2, 3, 4, 5, 55, 64, 65, 66, 67, 60 }; // zygomatic, infraorbital, parotid
    std::vector<int> leftZygoInfraParoRegionLandmarkIds{ 27, 28, 29, 30, 31, 32, 72, 73, 74, 75, 68, 59, 82 }; // zygomatic, infraorbital, parotid

    std::map<std::string, std::vector<int>> faceSkinSubRegions;


    faceSkinSubRegions["Nasal"] = nasalRegionLandmarkIds;
    faceSkinSubRegions["RightOrbital"] = rightOrbitalRegionLandmarkIds;
    faceSkinSubRegions["LeftOrbital"] = leftOrbitalRegionLandmarkIds;
    faceSkinSubRegions["Mental"] = mentalRegionLandmarkIds;
    faceSkinSubRegions["RightBuccal"] = rightBuccalRegionLandmarkIds;
    faceSkinSubRegions["LeftBuccal"] = leftBuccalRegionLandmarkIds;
    faceSkinSubRegions["RightZygoInfraParo"] = rightZygoInfraParoRegionLandmarkIds;
    faceSkinSubRegions["LeftZygoInfraParo"] = leftZygoInfraParoRegionLandmarkIds;

    std::map<std::string, int> faceSkinSubRegionsClassNumbers;

    faceSkinSubRegionsClassNumbers["Nasal"] = 20;
    faceSkinSubRegionsClassNumbers["RightOrbital"] = 21;
    faceSkinSubRegionsClassNumbers["LeftOrbital"] = 22;
    faceSkinSubRegionsClassNumbers["Mental"] = 23;
    faceSkinSubRegionsClassNumbers["RightBuccal"] = 24;
    faceSkinSubRegionsClassNumbers["LeftBuccal"] = 25;
    faceSkinSubRegionsClassNumbers["RightZygoInfraParo"] = 26;
    faceSkinSubRegionsClassNumbers["LeftZygoInfraParo"] = 27;




    std::map<std::string, cv::Mat> regionMasks;
    for (const auto& [regionName, landmarkIds] : faceSkinSubRegions)
    {
        OFIQ::Landmarks subRegionLandmarks;
        for (int i = 0; i < landmarkIds.size(); i++) {
            subRegionLandmarks.emplace_back(scaledLandmarks[landmarkIds[i]]);
        }

        cv::Mat mask = getPolygonMask(subRegionLandmarks, 200, 1);
        regionMasks.insert(make_pair(regionName, mask));
    }



    for (int i = 0; i < 200; i++) // For each pixel: 200 is the width and height of the image
    {
        for (int j = 0; j < 200; j++) // For each pixel: 200 is the width and height of the image
        {
            if (faceParsingImage.at<uchar>(i, j) == 1) // face skin is value 1
            {

                bool isPixelInASubRegion = false;
                for (const auto& [regionName, mask] : regionMasks)
                {
                    if (mask.at<uchar>(i, j) == 1) // we have chosen 1 as the pixel value for masks for now
                    {
                        faceParsingImage.at<uchar>(i, j) = faceSkinSubRegionsClassNumbers[regionName];
                        isPixelInASubRegion = true;
                        break; // We have the most important regions listed first, so the first one we encounter is the one that the pixel should belong to
                    }
                }

                if (!isPixelInASubRegion)
                {
                    faceParsingImage.at<uchar>(i, j) = 255; // residual face skin
                }

            }
        }
    }


    std::string inputFileWithoutExtension = inputFile.substr(0, inputFile.find_last_of("."));
    std::string outputPathFaceParsing = "./FGFP_images/FGFP_" + inputFileWithoutExtension + ".png"; // save as png because face parsing image needs to be lossless, since  the pixel value is the face parsing class
    cv::imwrite(outputPathFaceParsing, faceParsingImage);

    std::string outputPathAlignedImage = "./aligned_images/aligned_" + inputFile;
    cv::Mat alignedImage = session.getAlignedFace().clone();
    cv::resize(alignedImage, alignedImage, cv::Size(222, 222)); // face parsing image was scaled from (616-60)px to 200px. The aligned image must be scaled to 200+(200/(616-60))*60=222 (approx) to keep same ratio between face parsing image and aligned image.
    cv::imwrite(outputPathAlignedImage, alignedImage);
}


ReturnStatus OFIQImpl::vectorQuality(
    const OFIQ::Image& image, OFIQ::FaceImageQualityAssessment& assessments, const std::string& inputFile)
{
    auto session = Session(image, assessments);

    try
    {
        log("perform preprocessing:\n");
        performPreprocessing(session, inputFile);
    }
    catch (const OFIQError& e)
    {
        log("OFIQError: " + std::string(e.what()) + "\n");
        for (const auto& measure : m_executorPtr->GetMeasures() )
        {
            auto qualityMeasure = measure->GetQualityMeasure();
            switch (qualityMeasure)
            {
            case QualityMeasure::Luminance:
                session.assessment().qAssessments[QualityMeasure::LuminanceMean] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                session.assessment().qAssessments[QualityMeasure::LuminanceVariance] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                break;
            case QualityMeasure::CropOfTheFaceImage:
                session.assessment().qAssessments[QualityMeasure::LeftwardCropOfTheFaceImage] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                session.assessment().qAssessments[QualityMeasure::RightwardCropOfTheFaceImage] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                session.assessment().qAssessments[QualityMeasure::MarginBelowOfTheFaceImage] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                session.assessment().qAssessments[QualityMeasure::MarginAboveOfTheFaceImage] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                break;
            case QualityMeasure::HeadPose:
                session.assessment().qAssessments[QualityMeasure::HeadPoseYaw] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                session.assessment().qAssessments[QualityMeasure::HeadPosePitch] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                session.assessment().qAssessments[QualityMeasure::HeadPoseRoll] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                break;
            default:
                session.assessment().qAssessments[measure->GetQualityMeasure()] =
                { 0, -1, OFIQ::QualityMeasureReturnCode::FailureToAssess };
                break;
            }
        }

        return { e.whatCode(), e.what() };
    }

    log("execute assessments:\n");
    m_executorPtr->ExecuteAll(session);

    return ReturnStatus(ReturnCode::Success);
}

OFIQ_EXPORT std::shared_ptr<Interface> Interface::getImplementation()
{
    return std::make_shared<OFIQImpl>();
}

OFIQ_EXPORT void Interface::getVersion(int& major, int& minor, int& patch) const
{
    major = int(OFIQ_VERSION_MAJOR);
    minor = int(OFIQ_VERSION_MINOR);
    patch = int(OFIQ_VERSION_PATCH);
    return;
}
