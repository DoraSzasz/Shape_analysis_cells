#include "itkImage.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelToRGBImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkCustomColormapFunction.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include <sstream>
#include <cstdio>
#include <ctime>
#include <iostream>

const unsigned int Dimension = 3;
typedef unsigned char PixelType;
typedef unsigned long ConnCompType;
typedef float ParameterPixelType;
typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType, Dimension> ImageType;
typedef itk::Image<ParameterPixelType, Dimension> ParameterImageType;

typedef itk::Image<RGBPixelType, Dimension> RGBImageType;
typedef itk::Function::CustomColormapFunction<ImageType::PixelType, RGBImageType::PixelType> ColormapType;

typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;

typedef itk::Image<ConnCompType, Dimension> CCImageType;
typedef itk::ShapeLabelObject<ConnCompType,Dimension> ShapeLabelObjectType;
typedef itk::LabelMap<ShapeLabelObjectType> LabelMapType;

int main(int argc, char* argv[]) {

	std::clock_t start;
	double duration;

    	start = std::clock();
	
	PixelType distanceThreshold = 4;
	ImageType::Pointer image;
	if( argc < 2 ) {
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " inputImageFile" << std::endl;
		return EXIT_FAILURE;
	}
	else {
		typedef itk::ImageFileReader<ImageType> ReaderType;
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName(argv[1]);
		reader->Update();
		image = reader->GetOutput();
	}
	
	typedef itk::Image< unsigned char, Dimension > OutputImageType;
	typedef itk::ConnectedComponentImageFilter< ImageType, CCImageType > ConnectedComponentImageFilterType;
	typedef	itk::LabelImageToShapeLabelMapFilter< CCImageType, LabelMapType> I2LType;	
	typedef itk::RelabelComponentImageFilter<CCImageType, CCImageType> RelabelFilterType;
	
	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
	connected->SetInput(image);
	connected->SetFullyConnected(false);
	connected->Update();
	std::cout << "Number of objects: " << connected->GetObjectCount() << std::endl;
	
	RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
	relabelFilter->SetInput(connected->GetOutput());
	relabelFilter->Update();


	typedef std::vector< itk::SizeValueType > SizesInPixelsType ;
  const SizesInPixelsType & sizesInPixels = relabelFilter->GetSizeOfObjectsInPixels();
  SizesInPixelsType::const_iterator sizeItr = sizesInPixels.begin();
  SizesInPixelsType::const_iterator sizeEnd = sizesInPixels.end();
  std::cout << "Number of pixels per class " << std::endl;
  unsigned int kclass = 0;
  while (sizeItr != sizeEnd)
    {
    std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
    ++kclass;
    ++sizeItr;
    }

	I2LType::Pointer i2l = I2LType::New();
	i2l->SetInput( connected->GetOutput() );
	i2l->Update();


	LabelMapType *labelMap = i2l->GetOutput();
  // Retrieve all attributes
	std::ofstream outfile ("cell_attributes.csv");
	std::vector<ParameterPixelType> parameterTable;
	parameterTable.resize(labelMap->GetNumberOfLabelObjects());
  for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
    outfile << itk::NumericTraits<LabelMapType::LabelType>::PrintType(labelObject->GetLabel()) << ",";
    outfile << labelObject->GetNumberOfPixels() << ",";
    outfile << labelObject->GetPhysicalSize() << ",";
    outfile << labelObject->GetCentroid()[0] << ",";
    outfile << labelObject->GetCentroid()[1] << ",";
    outfile << labelObject->GetCentroid()[2] << ",";
    outfile << labelObject->GetPrincipalMoments()[0] << ",";
    outfile << labelObject->GetPrincipalMoments()[1] << ",";
    outfile << labelObject->GetPrincipalMoments()[2] << ",";
    outfile << labelObject->GetPrincipalAxes()[0][0] << ",";
    outfile << labelObject->GetPrincipalAxes()[0][1] << ",";
    outfile << labelObject->GetPrincipalAxes()[0][2] << ",";
    outfile << labelObject->GetPrincipalAxes()[1][0] << ",";
    outfile << labelObject->GetPrincipalAxes()[1][1] << ",";
    outfile << labelObject->GetPrincipalAxes()[1][2] << ",";
    outfile << labelObject->GetPrincipalAxes()[2][0] << ",";
    outfile << labelObject->GetPrincipalAxes()[2][1] << ",";
    outfile << labelObject->GetPrincipalAxes()[2][2] << ",";
    outfile << labelObject->GetElongation() << ",";
    outfile << labelObject->GetRoundness() << ",";
	parameterTable[n] = labelObject->GetNumberOfPixels();
    }

	CCImageType::RegionType inputRegion = connected->GetOutput()->GetRequestedRegion();
	ParameterImageType::Pointer parameterImage = ParameterImageType::New();
	parameterImage->SetRegions(inputRegion);

	parameterImage->SetSpacing(connected->GetOutput()->GetSpacing());
	parameterImage->SetOrigin(connected->GetOutput()->GetOrigin());
	parameterImage->Allocate();

  typedef itk::ImageRegionConstIterator< CCImageType > ConstIteratorType ;
  typedef itk::ImageRegionIterator< ParameterImageType > IteratorType; 

	ConstIteratorType inputIt(   connected->GetOutput(), inputRegion  );
	IteratorType      outputIt(  parameterImage,         inputRegion );
	inputIt.GoToBegin();
	outputIt.GoToBegin();
	while( !inputIt.IsAtEnd() )
		{
		float value= parameterTable[ inputIt.Get()];
		outputIt.Set(value==1? 255: 0); 
		++inputIt;
		++outputIt;
		}


	typedef itk::ImageFileWriter<ParameterImageType> WriterType;
        WriterType::Pointer writer2 = WriterType::New();
        writer2->SetFileName(argv[2]);
        writer2->SetInput(parameterImage);
        writer2->Update();
	

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    	std::cout<<"printf: "<< duration <<'\n';

	return EXIT_SUCCESS;
}
