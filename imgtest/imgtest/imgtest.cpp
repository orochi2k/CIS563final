// imgtest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <windows.h> // This is a must
#include "iostream"

using namespace std;


//RGBTRIPLE get_pixel(int x,int y)
//{
// // Image define from earlier
// return image[(bih.biHeight-1-y)*bih.biWidth+x];
//}
//RGBTRIPLE set_pixel(int x,int y,RGBTRIPLE color)
//{
// // Image define from earlier
// image[(bih.biHeight-1-y)*bih.biWidth+x] = color;
//}



int _tmain(int argc, _TCHAR* argv[])
{
	HANDLE hfile;
DWORD written;
BITMAPFILEHEADER bfh;
BITMAPINFOHEADER bih;
RGBTRIPLE *image;
hfile = CreateFile(L"C:/Users/orochi2k/CIS563final/test1.bmp",GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,NULL,NULL);
 // Read the header
 ReadFile(hfile,&bfh,sizeof(bfh),&written,NULL);
 ReadFile(hfile,&bih,sizeof(bih),&written,NULL);
 // Read image
 int imagesize = bih.biWidth*bih.biHeight; // Helps you allocate memory for the image
 image = new RGBTRIPLE[imagesize]; // Create a new image (I'm creating an array during runtime)
 ReadFile(hfile,image,imagesize*sizeof(RGBTRIPLE),&written,NULL); // Reads it off the disk
 // Close source file
 CloseHandle(hfile);
 // Now for some information
 cout<<"The image width is "<<bih.biWidth<<"\n"; // Will output the width of the bitmap
 cout<<"The image height is "<<bih.biHeight<<"\n"; // Will output the height of the bitmap
 RGBTRIPLE test;
 int xp,yp;
 xp = 16;
 yp = 4;
 test = image[(bih.biHeight-1-yp)*bih.biWidth+xp];
 cout<<"R:"<<(int)test.rgbtRed<<endl;
 cout<<"G:"<<(int)test.rgbtGreen<<endl;
 cout<<"B:"<<(int)test.rgbtBlue<<endl;
	/*Bitmap *image; 
	image=new Bitmap();
	if (image->loadBMP("C:/Users/orochi2k/CIS563final/test1.bmp")) {
		cout<<"w:"<<image->width<<" h:"<<image->height<<endl;
		
	}
	int w = 1;
	int h = 3;
	cout<<(int)image->data[(w + h * image->width) * 3 ]<<endl;
	cout<<(int)image->data[(w + h * image->width) * 3 + 1]<<endl;
	cout<<(int)image->data[(w + h * image->width) * 3 + 2]<<endl;*/



	system("pause");
	return 0;
}

