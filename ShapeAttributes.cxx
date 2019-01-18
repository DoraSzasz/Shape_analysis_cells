#include "itkImage.h"
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
typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType, Dimension> ImageType;

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
	
	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
	connected->SetInput(image);
	connected->SetFullyConnected(false);
	connected->Update();
	std::cout << "Number of objects: " << connected->GetObjectCount() << std::endl;

	I2LType::Pointer i2l = I2LType::New();
	i2l->SetInput( connected->GetOutput() );
	i2l->Update();


	LabelMapType *labelMap = i2l->GetOutput();
  // Retrieve all attributes
	std::ofstream outfile ("cell_attributes.txt");
  for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
    {
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
    outfile << "Label: "
              << itk::NumericTraits<LabelMapType::LabelType>::PrintType(labelObject->GetLabel()) << std::endl;
    outfile << "    BoundingBox: "
              << labelObject->GetBoundingBox() << std::endl;
    outfile << "    NumberOfPixels: "
              << labelObject->GetNumberOfPixels() << std::endl;
    outfile << "    PhysicalSize: "
              << labelObject->GetPhysicalSize() << std::endl;
    outfile << "    Centroid: "
              << labelObject->GetCentroid() << std::endl;
    outfile << "    NumberOfPixelsOnBorder: "
              << labelObject->GetNumberOfPixelsOnBorder() << std::endl;
    outfile << "    PerimeterOnBorder: "
              << labelObject->GetPerimeterOnBorder() << std::endl;
    outfile << "    FeretDiameter: "
              << labelObject->GetFeretDiameter() << std::endl;
    outfile << "    PrincipalMoments: "
              << labelObject->GetPrincipalMoments() << std::endl;
    outfile << "    PrincipalAxes: "
              << labelObject->GetPrincipalAxes() << std::endl;
    outfile << "    Elongation: "
              << labelObject->GetElongation() << std::endl;
    outfile << "    Perimeter: "
              << labelObject->GetPerimeter() << std::endl;
    outfile << "    Roundness: "
              << labelObject->GetRoundness() << std::endl;
    outfile << "    EquivalentSphericalRadius: "
              << labelObject->GetEquivalentSphericalRadius() << std::endl;
    outfile << "    EquivalentSphericalPerimeter: "
              << labelObject->GetEquivalentSphericalPerimeter() << std::endl;
    outfile << "    EquivalentEllipsoidDiameter: "
              << labelObject->GetEquivalentEllipsoidDiameter() << std::endl;
    outfile << "    Flatness: "
              << labelObject->GetFlatness() << std::endl;
    outfile << "    PerimeterOnBorderRatio: "
              << labelObject->GetPerimeterOnBorderRatio() << std::endl;
    }

	int no_cells = connected->GetObjectCount();
	GeneratorType::Pointer generator = GeneratorType::New();
	generator->Initialize();
	
	int R_color[no_cells - 1];
	int G_color[no_cells - 1];
	int B_color[no_cells - 1];
	for(int i = 0; i < (int)no_cells; i++) {
		R_color[i] = (int)generator->GetUniformVariate(0,255);
		G_color[i] = (int)generator->GetUniformVariate(0,255);
		B_color[i] = (int)generator->GetUniformVariate(0,255);
	}


	typedef itk::LabelToRGBImageFilter<CCImageType, RGBImageType> RGBFilterType;
	RGBFilterType::Pointer rgbFilter = RGBFilterType::New();
	rgbFilter->SetInput( connected->GetOutput());
	
	for(int i = 0; i < (int)no_cells; i++) {
		rgbFilter->AddColor(R_color[i], G_color[i], B_color[i]);
	}

	std::cout << "Number of colors: " << rgbFilter->GetNumberOfColors() << std::endl;

	/*typedef itk::ImageFileWriter<RGBImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(rgbFilter->GetOutput());
	writer->Update();*/

	typedef itk::ImageFileWriter<RGBImageType> WriterType;
        WriterType::Pointer writer2 = WriterType::New();
        writer2->SetFileName(argv[2]);
        writer2->SetInput(rgbFilter->GetOutput());
        writer2->Update();
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    	std::cout<<"printf: "<< duration <<'\n';

	return EXIT_SUCCESS;
}
