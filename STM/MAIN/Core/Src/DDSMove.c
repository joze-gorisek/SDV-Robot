#include "DDSM115.h"


float DDSangle01           = 0.0f;      //Current angle    [0-32767] -> [0 360]   address= 0x10 over RS485
float DDSrpm01             = 0.0f;      //Current velocity [0-     ] -> [-330RPM  330RPM]
float DDScurrent01         = 0.0f;      //Current Torque   [-32767   32767] -> [-8A   8A]
float DDSvelocityRadial01  = 0.0f;      //Velocity [rad/s]
float DDSvelocity01        = 0.0f;      // translational velocity [m/s]
float DDSangle01_old  = 0.0f;
float DDSangle01m_old = 0.0f;


float DDSangle30           = 0.0f;     // Current angle    [0-32767] -> [0 360]   address= 0x30 over RS485
float DDSrpm30             = 0.0f;     // Current velocity [0-     ] -> [-330RPM  330RPM]
float DDScurrent30         = 0.0f;     // Current Torque   [-32767   32767] -> [-8A   8A]
float DDSvelocityRadial30  = 0.0f;     // Velocity   [rad/s]
float DDSvelocity30        = 0.0f;     // Translational velocity [m/s]
float DDSangle30_old = 0.0f;
float DDSangle30m_old = 0.0f;

float Rwheel = 0.05; //Wheel Radius [m]
float delatx01=0.0;  //Initial position offset
float delatx30=0.0;  //Initial position offset


float distnacex01=0.0,anglex01=0.0,anglex01_old=0.0, anglex01_offset=0.0;
long int revolutionx01=0;//Number of wheel turns
int   zerocrossx01up=0;
int   zerocrossx01down=0;


float distnacex30=0.0,anglex30=0.0,anglex30_old=0.0, anglex30_offset=0.0;
long int revolutionx30=0;//Number of wheel turns
int   zerocrossx30up=0;
int   zerocrossx30down=0;

/*Measure distance of the robot from DDSM motor encoders*/
void Distnacex01()
{


	if( (DDSangle01m_old <= 360) && (DDSangle01m_old > 300) && (DDSangle01>=0) && (DDSangle01<100) )
		{
		  if(zerocrossx01up==0)
		  {
			  revolutionx01++;
		      zerocrossx01up=1;
		      zerocrossx01down=0;
		  }


		   //distnacex01 =2 * PI * Rwheel * revolutionx01  + (DDSangle01 * Rwheel);
		}
	else if( (DDSangle01m_old>=0) && (DDSangle01m_old<100) &&  (DDSangle01<=360) && (DDSangle01>300))
	   {
		      if(zerocrossx01down==0)
				  {
					  revolutionx01--;
				      zerocrossx01down=1;
				      zerocrossx01up=0;
				  }

			//revolutionx01--;
			//distnacex01 = 2 * PI * Rwheel * revolutionx01 - (360 - DDSangle01)* Rwheel;
	   }else if( (DDSangle01>100) && (DDSangle01_old<300))
	   {

		   zerocrossx01down=0;
		   zerocrossx01up=0;
	   }

	/*translational distnace [m]*/
	distnacex01 = -1.0 * (2 * PI * Rwheel * revolutionx01 +  DDSangle01* Rwheel * PI/180.0 - delatx01);  // delatx01 - compensate offset from encoder first start
	DDSangle01m_old = DDSangle01;

	/*translational velocity [m/s]*/
	DDSvelocity01 = DDSvelocityRadial01 * Rwheel;

}


/*Measure distance of the robot from DDSM motor encoders*/
void Distnacex30()
{


	if( (DDSangle30m_old <= 360) && (DDSangle30m_old > 300) && (DDSangle30>=0) && (DDSangle30<100) )
		{
		  if(zerocrossx30up==0)
		  {
			  revolutionx30++;
		      zerocrossx30up=1;
		      zerocrossx30down=0;
		  }


		   //distnacex01 =2 * PI * Rwheel * revolutionx01  + (DDSangle01 * Rwheel);
		}
	else if( (DDSangle30m_old>=0) && (DDSangle30m_old<100) &&  (DDSangle30<=360) && (DDSangle30>300))
	   {
		      if(zerocrossx30down==0)
				  {
					  revolutionx30--;
				      zerocrossx30down=1;
				      zerocrossx30up=0;
				  }

			//revolutionx01--;
			//distnacex01 = 2 * PI * Rwheel * revolutionx01 - (360 - DDSangle01)* Rwheel;
	   }else if( (DDSangle30>100) && (DDSangle30_old<300))
	   {

		   zerocrossx30down=0;
		   zerocrossx30up=0;
	   }

	/*translational distnace [m]*/
	distnacex30 = (2 * PI * Rwheel * revolutionx30 +  DDSangle30* Rwheel * PI/180.0 - delatx30);  //delatx01 - compensate offset from encoder first start
	DDSangle30m_old = DDSangle30;


	/*translational velocity [m/s]*/
    DDSvelocity30 = -1.0 * DDSvelocityRadial30 * Rwheel;

}
