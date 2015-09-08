/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include <iostream>
#include "itkImage.h"
#include "btkImageHelper.h"

float computeNLM(itk::Image< float, 3 >::Pointer & inputImage,int i, int j, int k, int hps, int hss, float h)
{
    itk::Image< float, 3 >::SizeType    size    = inputImage->GetLargestPossibleRegion().GetSize();

    itk::Image< float, 3 >::IndexType p;
    p[0] = i;
    p[1] = j;
    p[2] = k;

    int xmin = i - hss;
    if(xmin < hps) xmin = hps;
    int ymin = j - hss;
    if(ymin < hps) ymin = hps;
    int zmin = k - hss;
    if(zmin < hps) zmin = hps;
    int xmax = i + hss;
    if(xmax > size[0]-hps) xmax = size[0]-hps;
    int ymax = j + hss;
    if(ymax > size[1]-hps) ymax = size[0]-hps;
    int zmax = k + hss;
    if(zmax > size[2]-hps) zmax = size[0]-hps;

    float sum = 0;
    float outputValue = 0;
    for(int ii=xmin; ii<xmax; ii++)
        for(int jj=ymin; jj<ymax; jj++)
            for(int kk=zmin; kk<zmax; kk++)
            {
                float dist = 0;
                itk::Image< float, 3 >::IndexType pp;
                pp[0] = ii;
                pp[1] = jj;
                pp[2] = kk;

                for(int iii=-hps; iii< hps+1; iii++)
                    for(int jjj=-hps; jjj<hps+1; jjj++)
                        for(int kkk=-hps; kkk<hps+1; kkk++)
                        {
                            itk::Image< float, 3 >::IndexType p1;
                            p1[0] = p[0]+iii;
                            p1[1] = p[1]+jjj;
                            p1[2] = p[2]+kkk;
                            itk::Image< float, 3 >::IndexType p2;
                            p2[0] = pp[0]+iii;
                            p2[1] = pp[1]+jjj;
                            p2[2] = pp[2]+kkk;

                            dist += (inputImage->GetPixel(p1)-inputImage->GetPixel(p2))*(inputImage->GetPixel(p1)-inputImage->GetPixel(p2));
                        }
                float weight = exp (-dist/h);
                sum+= weight;
                outputValue+= weight*inputImage->GetPixel(pp);
            }
    if (sum>0)
        outputValue /= sum;
    return outputValue;
}
float computeNLMkji(itk::Image< float, 3 >::Pointer & inputImage,int i, int j, int k, int hps, int hss, float h)
{
    itk::Image< float, 3 >::SizeType    size    = inputImage->GetLargestPossibleRegion().GetSize();

    itk::Image< float, 3 >::IndexType p;
    p[0] = i;
    p[1] = j;
    p[2] = k;

    int xmin = i - hss;
    if(xmin < hps) xmin = hps;
    int ymin = j - hss;
    if(ymin < hps) ymin = hps;
    int zmin = k - hss;
    if(zmin < hps) zmin = hps;
    int xmax = i + hss;
    if(xmax > size[0]-hps) xmax = size[0]-hps;
    int ymax = j + hss;
    if(ymax > size[1]-hps) ymax = size[0]-hps;
    int zmax = k + hss;
    if(zmax > size[2]-hps) zmax = size[0]-hps;

    float sum = 0;
    float outputValue = 0;
    for(int kk=zmin; kk<zmax; kk++)
        for(int jj=ymin; jj<ymax; jj++)
            for(int ii=xmin; ii<xmax; ii++)
            {
                float dist = 0;
                itk::Image< float, 3 >::IndexType pp;
                pp[0] = ii;
                pp[1] = jj;
                pp[2] = kk;

                for(int kkk=-hps; kkk<hps+1; kkk++)
                    for(int jjj=-hps; jjj<hps+1; jjj++)
                        for(int iii=-hps; iii< hps+1; iii++)
                        {
                            itk::Image< float, 3 >::IndexType p1;
                            p1[0] = p[0]+iii;
                            p1[1] = p[1]+jjj;
                            p1[2] = p[2]+kkk;
                            itk::Image< float, 3 >::IndexType p2;
                            p2[0] = pp[0]+iii;
                            p2[1] = pp[1]+jjj;
                            p2[2] = pp[2]+kkk;

                            dist += (inputImage->GetPixel(p1)-inputImage->GetPixel(p2))*(inputImage->GetPixel(p1)-inputImage->GetPixel(p2));
                        }
                float weight = exp (-dist/h);
                sum+= weight;
                outputValue+= weight*inputImage->GetPixel(pp);
            }
    if (sum>0)
        outputValue /= sum;
    return outputValue;
}

int main( int argc, char * argv [] )
{
    std::string filename = "/Users/rousseau/Data/debug/newfbr/neo1.nii.gz";

    typedef itk::Image< float, 3 >  itkFloatImage;
    itkFloatImage::Pointer inputImage  = btk::ImageHelper< itkFloatImage > ::ReadImage(filename);
    itkFloatImage::Pointer outputImage = btk::ImageHelper< itkFloatImage > ::CreateNewImageFromPhysicalSpaceOf(inputImage,0.0);

    itkFloatImage::SizeType    size    = inputImage->GetLargestPossibleRegion().GetSize();

    int ps = 3; // patch size
    int hps= 1; //half patch size
    int hss= 2; //half search size

    float h = 10000;
    int i,j,k;
    //#pragma omp parallel for private(i,j,k) schedule(dynamic)
    /*
    for(i=hps; i<size[0]-hps; i++)
        for(j=hps; j<size[1]-hps; j++)
            for(k=hps; k<size[2]-hps; k++)
            {
                itk::Image< float, 3 >::IndexType p;
                p[0] = i;
                p[1] = j;
                p[2] = k;

                float outputValue = computeNLM(inputImage, i,j,k,hps,hss,h);
                outputImage->SetPixel(p,outputValue);
            }
    */
    #pragma omp parallel for private(i,j,k) schedule(dynamic)
    for(k=hps; k<size[2]-hps; k++)
        for(j=hps; j<size[1]-hps; j++)
            for(i=hps; i<size[0]-hps; i++)
            {
                itk::Image< float, 3 >::IndexType p;
                p[0] = i;
                p[1] = j;
                p[2] = k;

                float outputValue = computeNLMkji(inputImage, i,j,k,hps,hss,h);
                outputImage->SetPixel(p,outputValue);
            }


    btk::ImageHelper<itkFloatImage>::WriteImage(outputImage, "btktoto.nii.gz");



}
