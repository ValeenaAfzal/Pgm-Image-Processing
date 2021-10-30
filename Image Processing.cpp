#include<iostream>
#include<fstream>
#include<string>
#include <cmath>
#include<cstdlib>
#include<math.h>
using namespace std;
#define MaxRows 500
#define MaxCols 500
int count=0,count1=0;
// Code By Valeena Afzal 20L-1035 1B1 and Rameen Amir 20L-1283 1B1
struct grayImage
{
	grayImage()
    {
        Rows = Cols = 0;
        Loaded = false;
        for(int r = 0; r< MaxRows; r++)
            for(int c = 0; c< MaxCols; c++)
                Image[r][c] = 0;
    }
    unsigned short setPixel(unsigned short value, int r, int c)
    {
        if( r >= Rows || c >= Cols || r < 0 || c < 0)
            return -1;
        Image[r][c] = value;
        return value;
    }
    int getPixel(int r, int c)
    {
        if( r >= Rows || c >= Cols || r < 0 || c < 0)
            return -1;
        return Image[r][c];
    }
    int setRows(int rows)
    {
        if(rows < 1 || rows > MaxRows)
            return -1;
        Rows = rows;
        return Rows;
    }
    int getRows()
    {
        return Rows;
    }
    int setCols(int cols)
    {
        if(cols < 1 || cols > MaxCols)
            return -1;
        Cols = cols;
        return Cols;
    }
    int getCols()
    {
        return Cols;
    }
    int Save(string File_Name)
	{
    	int n=File_Name.length();
    	string e=File_Name.substr(n-4,4);
    	if (e.compare(".pgm") !=0 )
    	{
    		cout<<"File format is not correct"<<endl;
    		return 1;
		}
			ofstream OUT (File_Name.c_str());
			{
				if(OUT)
				{
					cout<<"Saved Successfully"<<endl;
			    	OUT<< "P2\n";
			    	OUT<< "# The optional comment by Valeena and Rameen\n ";
					OUT<<Cols<<" "<<Rows<<endl<<Maximum<<endl;
			    	for (int r=0;r<Rows;r++)
			    	{
			    	    for (int c=0;c<Cols;c++)
			    	    OUT<<Image[r][c]<<" ";
			    	    OUT<<endl;
			    	}
			    	
			    	OUT.close();
			    	return 0;
				}
			}
			if(!OUT)
			{
				cout<<"Unable to Create a file"<<endl;
				return 2;
			}
	}
    int load(string File_Name)
	{
    	int n=File_Name.length();
    	string e=File_Name.substr(n-4,4);
    	if (e.compare(".pgm") !=0 )
    	{
    		cout<<"File format is not correct";
    		return 1;
		}
    	char p[200];
    	char s[10];
    	ifstream IN(File_Name.c_str());
    	if(IN)
    	{
		cout<<"Loaded Successfully"<<endl;
		IN.getline(s,10);
    	IN.getline(p,200);
    	int rows,cols,maxval;
    	IN>>cols>>rows>>maxval;
    	setRows(rows);
		setCols(cols);
    	Maximum=maxval;
		Rows=rows;
		Cols=cols;
    	for (int i=0;i<rows;i++)
		{
			for (int j=0; j<cols; j++)
			{
				int val;
				IN>>val;
				Image[i][j]=val;
			}
		}
			Loaded=true;
		    return 0;
		}
		else
		{
			cout<<"Not Loaded"<<endl;
			return 2;
		}
	}
    void combineSideBySide(grayImage &Two, int fillValue = 0)
	{
		if (Rows>Two.Rows && Cols>Two.Cols)
		{
			count1=1;
			int L=0,T=0,R=0,B=0;
	    	for (int i=0;i<Two.Rows;i++)
	    	{
	    		for (int j=0;j<Two.Cols;j++)
	    		{
	    			if(Cols+j<MaxCols)
					Image[i][Cols+j]=Two.Image[i][j];
				}
	
			}
			T=Two.Rows;
			B=Rows;
			L=Cols;
			R=MaxCols;
			Fill(L,T,R,B, fillValue);
			if (Rows<Two.Rows)
			Rows=Two.Rows;
			Cols+=Two.Cols;
			if(Cols>MaxCols)
				Cols=MaxCols;
		}
		if (Two.Rows>Rows && Two.Cols>Cols) 
		{
			count1=2;
			int L=0,T=0,R=0,B=0;
	    	for (int i=0;i<Rows;i++)
	    	{
	    		for (int j=0;j<Cols;j++)
	    		{
	    			if(Two.Cols+j<MaxCols)
					Two.Image[i][Two.Cols+j]=Image[i][j];
				}
	
			}
			T=Rows;
			B=Two.Rows;
			L=Two.Cols;
			R=MaxCols;
			for(int i = T; i<= B; i++)
            {
            	for(int j = L; j <= R; j++)
                	Two.Image[i][j] = fillValue;
			}
			if (Two.Rows<Rows)
			Two.Rows=Rows;
			Two.Cols+=Cols;
			if(Two.Cols>MaxCols)
				Two.Cols=MaxCols;
		}

	}
    void combineTopToBottom(grayImage &Two, int fillValue = 0)
	{
		if (Rows>Two.Rows && Cols>Two.Cols)
		{
			count=1;
			int L=0,T=0,R=0,B=0;
			for (int i=0;i<Two.Rows;i++)
			{
				for (int j=0;j<Two.Cols;j++)
				{
					if(Rows+i<MaxRows)
					Image[Rows+i][j]=Two.Image[i][j];
				}
	
			}
			T=Rows;
			B=MaxRows;
			L=Two.Cols;
			R=Cols;
			Fill(L,T,R,B, fillValue);
			if (Cols<Two.Cols)
				Cols=Two.Cols;
				Rows+=Two.Rows;
			if(Rows>MaxRows)
				Rows=MaxRows;
		}
		if (Two.Rows>Rows && Two.Cols>Cols)
		{
			count=2;
			int L=0,T=0,R=0,B=0;
			for (int i=0;i<Rows;i++)
			{
				for (int j=0;j<Cols;j++)
				{
					if(Two.Rows+i<MaxRows)
					Two.Image[Two.Rows+i][j]=Image[i][j];
				}
	
			}
			T=Two.Rows;
			B=Rows+Two.Rows;
			L=Cols;
			R=Two.Cols;
			for(int i = T; i<= B; i++)
            {
            	for(int j = L; j <= R; j++)
                	Two.Image[i][j] = fillValue;
			}
			
			if (Cols>Two.Cols)
				Two.Cols=Cols;
				Two.Rows+=Rows;
			if(Rows<MaxRows)
				Two.Rows=MaxRows;
		}

    }
    void FadeIn (grayImage &second,int seconds, int frameRate , string baseFileName)
	{
		
		double size=1.0/(frameRate*seconds);
		int ro=Rows,co=Cols;
		if (ro<second.Rows)
		ro=second.Rows;
		if (co<second.Cols)
		co=second.Cols;
		Loaded=true;
		Maximum=Maximum;
		Rows=ro;
		Cols=co;
		if (Maximum<second.Maximum)
		Maximum=second.Maximum;
		int counter=0;
		for (double a=1;a>=0 ;a-=size)
		{
			for(int i=0;i<ro;i++)
			{
				for(int j=0;j<co;j++)
				{
					Image[i][j]=Image[i][j]*a+(1-a)*second.Image[i][j];

				}
			}
			int gg=10; char V[gg];
			itoa(counter, V , gg);
			Save(baseFileName+V+".pgm");
			counter++;
			if (0>a-size && a-size> -size)
				a=0;

		}
	}
    void Rotate(grayImage& RotatedImage, double angle = 90, int aboutx = 0, int abouty = 0 )
	{
		int rr,cc;
		angle=(3.141592*angle)/180;
		for(int i=0;i<Rows;i++)
		{
			for(int j=0;j<Cols;j++)
			{
				rr=(i-aboutx)*cos(angle) - (j-abouty)*sin(angle);
				cc=(j-abouty)*cos(angle) + (i-aboutx)*sin(angle);
				if(rr>0 && rr<MaxRows && cc>0 && cc<MaxCols)
					RotatedImage.Image[rr][cc]=Image[i][j];
			}
		}
		
		RotatedImage.Rows=MaxRows;
		RotatedImage.Cols=MaxCols;
		RotatedImage.Maximum=Maximum;
		RotatedImage.Loaded=true;
		RotatedImage.meanFilter(RotatedImage,3);

    }
    void Flip(grayImage & Result, int flip=0)
	{
  		if (flip)
  		{
  			FlipVertical(Result);
		}
		else
		{
			 FilpHorizontal(Result);
		}
	}
    void Negative()
	{
		int max=255,min=0;
		int lar=0;
			for (int r=0;r<Rows;r++)
            {
                for (int c=0;c<Cols;c++)
                {
                    if (Image[r][c]>lar)
                    lar=Image[r][c];
                }
            }
		  	for (int r=0;r<Rows;r++)
            {
                for (int c=0;c<Cols;c++)
                {
                    Image[r][c]-=lar;
                    Image[r][c]*=(-1);
                    if (Image[r][c]>max)
                        Image[r][c]=max;
                    if (Image[r][c]<min)
                        Image[r][c]=min;

                }
            }
	}
    void changeBrightness(int amount)
	{
		int max=255,min=0;
		for (int r=0;r<Rows;r++)
	    {
	        for (int c=0;c<Cols;c++)
	        {
				Image[r][c]+=amount;
		        if (Image[r][c]>max)
                    Image[r][c]=max;
			    if (Image[r][c]<min)
                    Image[r][c]=min;
			}
	    }
	}
    void Quantize(grayImage& Result,int Levels)
	{
        int step = (Maximum+1)/Levels;
        int c=step;
        int startpixel=0,endpixel=0;

			while(c<=256)
			{
				startpixel = c - step;
		        endpixel = c - 1;
				for(int i=0;i<Rows;i++)
				{
					for(int j=0;j<Cols;j++)
					{
						int value = (Image[i][j] / step )*step + step/2;
						Result.Image[i][j]=value;
					}
				}
			c+=step;
			}
			Result.Cols=Cols;
			Result.Rows=Rows;
			Result.Maximum=Maximum;
			Result.Loaded=true;
    	return ;
	}
    void medianFilter(grayImage& Result, int filterSize = 3)
	{
		int filter=filterSize*filterSize,u=0,arr[filter]={};
        for(int i=0;i<Rows-1;i++)
        {
        	for(int j=0;j<Cols-1;j++)
        	{
        		for(int a=0;a<filterSize;a++)
				{
        			for(int b=0;b<filterSize;b++)
					{
						int g=i-filterSize/2,h=j-filterSize/2;
						arr[u]=Image[a+g][b+h];
						u++;
					}
				}
				selectionSort(arr,0,filter);
				int r=(filterSize*filterSize)/2;
				Result.Image[i][j]=arr[r];
				u=0;
				for (int empty=0; empty<filter; empty++)
				{
					arr[empty]=0;
				}
				
        	}
        }
        Result.Rows=Rows;
        Result.Cols=Cols;
        Result.Maximum=Maximum;
        Result.Loaded=true;
    }
    void meanFilter(grayImage& Result, int filterSize = 3)
	{
		int n=0,c=0,sum=0,avg=0,k=1;
        if(filterSize!=3) k=(filterSize/3)+1;
        c=filterSize;
        n=filterSize*filterSize;
		for(int i=1;i<Rows-1;i++)
		{

			for(int j=1;j<Cols-1;j++)
			{
			sum=0;

				for(int a=i-k;a<=i+k;a++)
				{
					for(int b=j-k;b<=j+k;b++)
						sum+=Image[a][b];
				}

				avg=sum/n;
				Result.Image[i][j]=avg;
				avg=0;
			}
		}
		Result.Rows=Rows;
        Result.Cols=Cols;
        Result.Maximum=Maximum;
        Result.Loaded=true;

    }
    void Resize(grayImage& Result,int NewRows, int NewColumns)
	{
		double rows=NewRows,cols=NewColumns,rr=Rows,cc=Cols;
		double RR;
		RR=rows/rr;
		double CR;
		CR=cols/cc;
		int I=0,J=0;
		for (int i=0;i<NewRows;i++)
		{
			for (int j=0;j<NewColumns;j++)
			{
				I=i / RR;
				J=j / CR;
				Result.Image[i][j]=Image[I][J];
			}
		}
		
		Result.Rows=NewRows;
		Result.Cols=NewColumns;
		Result.Maximum=Maximum;
		Result.Loaded=true;
    }
    void Resize(grayImage& Result,double Ratio)
	{
		int NewRows=Rows*Ratio,NewColumns=Cols*Ratio;
		Resize(Result,NewRows,NewColumns);
	}
    void Transform(grayImage& Result,double Matrix[3][3])
	{
		int a=0,b=1,d=2;
		for (int r=0;r<Rows;r++)
		{
			for(int c=0;c<Cols;c++)
			{
				int I,J,K;
				I=(r*Matrix[a][a])+(c*Matrix[a][b])+Matrix[a][d];
				J=(r*Matrix[b][a])+(c*Matrix[b][b])+Matrix[b][d];
				K=(r*Matrix[d][a])+(c*Matrix[d][b])+Matrix[d][d];
				if (K!=0)
				{
					int NI=I/K, NJ=J/K;
					if(NI>=0 && NI<MaxRows && NJ>=0 && NJ< MaxCols)
					{
						Result.Image[NI][NJ]=Image[r][c];
					}
				}
			}
		}
		Result.Cols=Cols;
		Result.Rows=Rows;
		Result.Maximum=Maximum;
		Result.Loaded=true;
    }
    void Filter(grayImage& Result,double Mask[3][3])
	{
		int S=0,size=3;
        for(int i=0;i<Rows;i++)
        {
        	for(int j=0;j<Cols;j++)
        	{
        		S=0;
        		int I=i- size/2;
        		int J=j- size/2;
        		for(int r=0;r<size;r++)
				{
        			for(int c=0;c<size;c++)
        				S+=(Image[r+I][c+J]*Mask[r][c]);
        		}
				Result.Image[i][j]=S;
        	}
        }
		Result.Cols=Cols;
		Result.Rows=Rows;
		Result.Maximum=Maximum;
		Result.Loaded=true;
	}
    void DerivativeImage(grayImage& Result)
	{
		double Mask[3][3]={-1,0,1,-1,0,1,-1,0,1};
		double Mask2[3][3]={-1,-1,-1,0,0,0,1,1,1};
		                
        for(int i=0;i<=Rows-3;i++)
        {
        	for(int j=0;j<=Cols-3;j++)
        	{
        	    int S1 = 0 , S2 = 0, size=3;
        		for(int r=0;r<size;r++)				
        			for(int c=0;c<size;c++)
					{
        				S1+=( Image[r+i][c+j] * Mask[r][c] );				
        				S2+=( Image[r+i][c+j] * Mask2[r][c] );       				
					}
				int x = ( (S1*S1) + (S2*S2) );
				if (x<0) x=x*(-1);				
				int y = sqrt(x);
				Result.Image[i][j]=y;
        	}
        }

        Result.Cols=Cols;
		Result.Rows=Rows;
		Result.Maximum=Maximum;
		Result.Loaded=true;
    }
	void Crop(grayImage& Result,int L,int R,int T,int B,int RSF=0)
	{
		double Ratio;
		int cols=R-L;
		int rows=B-T;
		for(int i=T,a=0;i<=B,a<rows;i++,a++)
			for(int j=L,b=0;j<=R,b<cols;j++,b++){
				int val = Image[i][j];
				Result.Image[a][b]=val;
			}
		if (RSF==0)
		{
			Result.Rows=rows;
	        Result.Cols=cols;
		}
		else if(RSF==1)
		{
			cout<<"Enter New rows For resize: ";cin>>rows;
			cout<<"Enter New columns for resize: ";cin>>cols;
			Resize(Result,rows,cols);
		}
		else if(RSF==2)
		{
			cout<<"Enter ratio For resize: ";cin>>Ratio;
			Resize(Result,Ratio);
		}
		Result.Maximum=Maximum;
	    Result.Loaded=true;
	}
private:
	void FilpHorizontal(grayImage&Copy)
	{
		for (int r=0,k=Rows-1;r<Rows,k>-1;r++,k--)
	    {
	        for (int c=0;c<Cols;c++)
	        Copy.Image[r][c]=Image[k][c];
	    }
	    Copy.Cols=Cols;
		Copy.Rows=Rows;
		Copy.Maximum=Maximum;
		Copy.Loaded=true;
	}
	void FlipVertical(grayImage & copy)
	{
		for (int r=0;r<Rows;r++)
	    {
	        for (int c=0, k=Cols-1;c<Cols,k>-1;c++,k--)
	        copy.Image[r][c]=Image[r][k];
	    }
		copy.Cols=Cols;
		copy.Rows=Rows;
		copy.Maximum=Maximum;
		copy.Loaded=true;
    }
	void Fill(int L, int T, int R, int B, int FillValue)
	{
        for(int i = T; i<= B; i++)
            for(int j = L; j <= R; j++)
                Image[i][j] = FillValue;
    }
	void swapping(int &a, int &b)
	{
	    int temp;
	    temp = a;
	    a = b;
	    b = temp;
	}
	void display(int arr[],int SI,int EI)
	{
	    for(int i = SI; i<=EI; i++)
	      cout << arr[i] << " ";
		cout << endl;
	}
	int selectionSort(int arr[],int SI,int EI)
	{
		int i,j,imin;
		for(int i=SI;i<EI;i++)
		{
			imin=i;
			for(int j=i+1;j<=EI;j++)
				if(arr[j]<arr[imin])
					imin=j;
			swapping(arr[i],arr[imin]);
		}

	}
	int search1(int r,int arr[])
	{
		int x;
		for(int i=0;i<=r;i++)
		{
			if (i==r)
				return x=arr[r],x;
		}

	}

    unsigned short Image[MaxRows][MaxCols];
    int Rows, Cols, Maximum;
    bool Loaded;
};
    int getChoice(string Message, int Low_Limit, int Up_Limit)
	{
    int Choice = Low_Limit;
				
	do
	{
		cout<<Message<<" "<<Low_Limit<<" and "<<Up_Limit<<" ";
        cin>>Choice;
    }
	while(Choice < Low_Limit || Choice >Up_Limit);

    return Choice;
	}


int main()
{
    cout<<"      Hello   \n"<<endl;
    int NewRows=0, NewCols=0, Choice=0;
    double Ratio=0;
    //creating variables of type struct
    grayImage GM, GM2, GM3, GM4;
    double Mask[3][3];
    double Matrix[3][3];
    cout<<"Enter [3x3] Matrix for Transformation for futhur use";
	    for (int i=0;i<3;i++)
	    	{
	    		for (int j=0;j<3;j++)
	    			cin>>Matrix[i][j];
			}
		cout<<"Enter [3x3] Mask For Filter to use futhur ";
		for (int i=0;i<3;i++)
	    	{
	    		for (int j=0;j<3;j++)
	    			cin>>Mask[i][j];
			}
	string Menu = "\n\n 0. Load File\n 1. Save File\n 2. Combine Side by Side\n 3. Combine Top to Bottom\n 4. Fade In\n 5. Rotate\n 6. Flip\n 7. Negative\n 8. Change Brightness\n 9. Quantize\n10. Mean Filter\n11. Median Filter\n12. Transform\n13. Filter\n14. Resize with Columns and Rows\n15. Resize With Ratio\n16. Derivative\n17. Crop\n18. exit";
	do
	{
		
		cout<<Menu;
		cout<<"\n\nPlease enter the option number to perform an operation: ";
		cin>>Choice;
		if(Choice<0 || Choice > 18)
		{
			cout<<"\n\n\nArray NOT AVAILABLE\n\n\n";
			cin>>Choice;
		} 
		if(Choice == 0)
		{
			string File_Name;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
		}
		else if(Choice == 1)
		{
			string File_Name;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Name;
			GM.Save(File_Name);
		}		
		else if(Choice == 2)
		{
			string File_Name,File_Names;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter Second File Name With .pgm Extension To Load:  ";
			cin>>File_Names;
			GM2.load(File_Names);
			int fill;
			cout<<"Enter Fill Value";
			cin>>fill;
			GM.combineSideBySide(GM2,fill);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			if (count1==1)
			GM.Save(File_Namez);
			if (count1==2)
			GM2.Save(File_Namez);
				
		}		
		else if(Choice == 3)
		{
			string File_Name,File_Names;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter Second File Name With .pgm Extension To Load:  ";
			cin>>File_Names;
			GM2.load(File_Names);
			int fill;
			cout<<"Enter Fill Value";
			cin>>fill;
			GM.combineTopToBottom(GM2,fill);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			if (count==1)
			GM.Save(File_Namez);
			if (count==2)
			GM2.Save(File_Namez);	
		}		
		else if(Choice == 4)
		{
			string File_Name,File_Names,Base;
			int sec,frame;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter Second File Name With .pgm Extension To Load: ";
			cin>>File_Names;
			GM2.load(File_Names);
			cout<<"Enter BaseFileName: ";
			cin>>Base;
			cout<<"Enter Number of Seconds Per Image ";
			cin>>sec;
			cout<<"Enter Number of Frame ";
			cin>>frame;
			GM.FadeIn(GM2,sec,frame,Base);
		}		
		else if(Choice == 5)
		{
			string File_Name;
			int x,y;
			double angle;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter x-axis: ";
			cin>>x;
			cout<<"Enter y-axis: ";
			cin>>y;
			cout<<"Enter Angle: ";
			cin>>angle;
			GM.Rotate(GM2,angle,x,y);
			GM2.meanFilter(GM4,3);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM4.Save(File_Namez);
		}		
		else if(Choice == 6)
		{
			string File_Name;
			int c; 
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter 0 To Flip Vertical And Non Zero Number To Horizontal:  ";
			cin>>c;
			GM.Flip(GM2,c);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM2.Save(File_Namez);
		}		
		else if(Choice == 7)
		{
			string File_Name;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			GM.Negative();
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM.Save(File_Namez);	
		}		
		else if(Choice == 8)
		{
			string File_Name;
			int c;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"\nEnter amount for Brightness : ";
			cin>>c;
			GM.changeBrightness(c);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM.Save(File_Namez);	
		}		
		else if(Choice == 9)
		{
			string File_Name;
			int c;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter Level Of Quantization";
			cin>>c;
			GM.Quantize(GM2,c);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM2.Save(File_Namez);
		}
		else if(Choice == 10)
		{
			string File_Name;
			int c;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter Filter Size ";
			cin>>c;
			GM.meanFilter(GM4,c);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM4.Save(File_Namez);	
		}
		else if(Choice == 11)
		{
			string File_Name;
			int c;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter Filter Size ";
			cin>>c;
			GM.medianFilter(GM4,c);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM4.Save(File_Namez);
		}
		else if(Choice == 12)
		{
			string File_Name;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			GM.Transform(GM4,Matrix);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM4.Save(File_Namez);
		}
		else if(Choice == 13)
		{
			string File_Name;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			GM.Filter(GM2,Mask);
	    	string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM2.Save(File_Namez);
		}
		else if(Choice == 14)
		{
			string File_Name;
			int NewRows,NewCols;
			cout<<"Enter File Name With .pgm Extension To Load ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter new rows for image: ";
			cin>>NewRows;
			cout<<"Enter new columns for image: ";
			cin>>NewCols;
			GM.Resize(GM2,NewRows,NewCols);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM2.Save(File_Namez);
		}
		else if(Choice == 15)
		{
			string File_Name;
			double Ratio;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			cout<<"Enter Ratio: ";
			cin>>Ratio;
			GM.Resize(GM2,Ratio);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM2.Save(File_Namez);
		}
		else if(Choice == 16)
		{
			string File_Name;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			GM.load(File_Name);
			GM.DerivativeImage(GM2);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM2.Save(File_Namez);
		}
		else if(Choice == 17)
		{
			int L,R,T,B,RSF;
			string File_Name;
			cout<<"Enter File Name With .pgm Extension To Load:  ";
			cin>>File_Name;
			cout<<"Enter value for Left: ";
			cin>>L;
			cout<<"Enter value for Right: ";
			cin>>R;
			cout<<"Enter value for Top: ";
			cin>>T;
			cout<<"Enter value for Bottom: ";
			cin>>B;
			cout<<"Enter 0 if you want cropped size image: \nEnter 1 if you want to resize Rows Columns: \nEnter 2 if you want to resize with ratio :\n Enter Your Choice: :  ";cin>>RSF;
			GM.load(File_Name);			
			GM.Crop(GM2,L,R,T,B,RSF);
			string File_Namez;
			cout<<"Enter File Name With .pgm Extension To Save:  ";
			cin>>File_Namez;
			GM2.Save(File_Namez);
		}
	}
	while(Choice!=18);
	cout<<"\nThanks ! Hope that you are satisfied";
	return 0;
}
