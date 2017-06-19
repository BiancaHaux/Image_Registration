#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkTranslationTransform.h"
#include "itkRescaleIntensityImageFilter.h"


//Metriken zur Auswahl:
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"

//Main

int main(int argc, char *argv[])
{

	typedef itk::Image<double, 2>																	DoubleImageType;
	typedef itk::Image<unsigned char, 2>															UnsignedCharImageType;
	typedef itk::ImageFileReader<UnsignedCharImageType>												ReaderType;
	typedef itk::RescaleIntensityImageFilter< DoubleImageType, UnsignedCharImageType >				RescaleFilterType;
	typedef itk::ImageFileWriter< UnsignedCharImageType >											WriterType;
	

	if (argc < 3)
	{
		std::cout << "Usage: " << argv[0] << "imageFile1 imageFile2" << std::endl;
		return EXIT_FAILURE;
	}


	//image: Bild, das die Metrikwerte (Typ double) aufnimmt
	DoubleImageType::Pointer image = DoubleImageType::New(); 

	itk::Index<2> start;
	start.Fill(0);
	itk::Size<2> size;
	size[0] = 201;
	size[1] = 201;
	itk::ImageRegion<2> region(start, size);
	image->SetRegions(region);
	image->Allocate();


	//Bilder einlesen
	ReaderType::Pointer fixedReader = ReaderType::New();
	fixedReader->SetFileName(argv[1]);
	fixedReader->Update();

	ReaderType::Pointer movingReader = ReaderType::New();
	movingReader->SetFileName(argv[2]);
	movingReader->Update();

	UnsignedCharImageType::Pointer fixedImage = fixedReader->GetOutput();
	UnsignedCharImageType::Pointer movingImage = movingReader->GetOutput();

	//Fall MeanSquares Metrik
	typedef itk::MeanSquaresImageToImageMetric < UnsignedCharImageType, UnsignedCharImageType > MetricType;

	//Fall NormalizedCorrelation Metrik
	//typedef itk::NormalizedCorrelationImageToImageMetric < UnsignedCharImageType, UnsignedCharImageType > MetricType;


	//Interpolation und Transformation festlegen
	typedef itk::NearestNeighborInterpolateImageFunction<UnsignedCharImageType, double >		InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetInputImage(fixedImage);

	typedef itk::TranslationTransform < double, 2 >													TransformType;
	TransformType::Pointer transform = TransformType::New();


	//Metrik instantiieren und Prameter setzen
	MetricType::Pointer metric = MetricType::New();
	metric->SetFixedImage(fixedImage);
	metric->SetMovingImage(movingImage);
	metric->SetFixedImageRegion(fixedImage->GetLargestPossibleRegion());
	metric->SetTransform(transform);
	metric->SetInterpolator(interpolator);

	TransformType::ParametersType params(transform->GetNumberOfParameters());
	params.Fill(0.0);

	metric->Initialize();
	itk::Index<2> pixel;
	pixel[0]  =  0;
	pixel[1]  =  0;
	float wert;

	for (double x = -100.0; x <= 100.0; x += 1.0)
	{
		params(0) = x;
		
		for (double y = -100.0; y <= 100.0; y += 1.0)
		{
			params(1) = y;
			
			wert = metric->GetValue(params);
			
			image->SetPixel(pixel, wert);

			pixel[1]++;
		}
		pixel[1] = 0;
		pixel[0]++;
	}
	

	//Nun die Grauwerte aus image1 auf einen Bereich von 0 bis 255 reskalieren
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(image);
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();

	
	
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("metric.png");
	writer->SetInput(rescaleFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}
