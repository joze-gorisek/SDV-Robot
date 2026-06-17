/*
 *
 *
 * */

#include "CyberGear.h"


float front_angle=     -0.0f;
float back_angle=      -0.0f;
float left_angle=      -0.0f;
float right_angle=     -0.0f;

float front_angle_old=  -0.0f;
float back_angle_old=   -0.0f;
float left_angle_old=   -0.0f;
float right_angle_old=  -0.0f;

float roll_right     =0.0f;
float roll_left      =0.0f;
float roll_right_old =0.0f;
float roll_left_old  =0.0f;



float Up_angle =       0.0f;
float Up_angle_old =   0.0f;

float speedLimit     = 3.5f; // e.g. 1 rad/s  max 30rad/s
float speedLimit_old = 0.8f;
int  cntflag=0;
extern int Cybercounter;

/*Main Cyber-movements function*/
void BodyMove()
{

	if(Up_angle-Up_angle_old!=0)
	{
		CyberUpDown();

	}else if(front_angle-front_angle_old!=0)
	{
		CyberFrontAngle();
		/*Update values*/
		front_angle_old=front_angle;
		Up_angle       =front_angle;
		Up_angle_old   =front_angle;
	    roll_left_old  =front_angle;
	    roll_right_old =front_angle;
	    roll_left      =front_angle;
	  	roll_right     =front_angle;

	}else if(back_angle-back_angle_old!=0)
	{

		CyberBackAngle();
		/*Update values*/
		back_angle_old =back_angle;
		Up_angle       =back_angle;
		Up_angle_old   =back_angle;
	    roll_left_old  =back_angle;
	    roll_right_old =back_angle;
	    roll_left      =back_angle;
	  	roll_right     =back_angle;

	}else if(speedLimit-speedLimit_old!=0)
	{
		SpeedLimit();

	}else if(roll_right-roll_right_old!=0.0)
	{

		CyberRollRight();
		roll_right_old  =roll_right;
		back_angle_old  =roll_right;
		back_angle      =roll_right;
		front_angle_old =roll_right;
	    front_angle     =roll_right;
	    Up_angle        =roll_right;
	    Up_angle_old    =roll_right;

	}else if(roll_left-roll_left_old!=0.0)
	{

		CyberRollLeft();
		    roll_left_old   =roll_left;
			back_angle_old  =roll_left;
			back_angle      =roll_left;
			front_angle_old =roll_left;
		    front_angle     =roll_left;
		    Up_angle        =roll_left;
		    Up_angle_old    =roll_left;

	}




}



/*Cyber front-angle settings*/
void CyberFrontAngle(void)
{
	//if(Cybercounter==0)     /*Check if Cybergear is ready for control*/
//	{
	      if(front_angle>1.7)
	    	  front_angle=1.7;
	      else if(front_angle<0)
	    	  front_angle=0.0;

		  SetAngle(front_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_1_ID);
		  SetAngle(-front_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_2_ID);
	//}

}

/*Cyber back-angle settings*/
void CyberBackAngle(void)
{
	  if(back_angle>1.7)
		  back_angle=1.7;
	  else if(back_angle<0)
		  back_angle=0.0;

		SetAngle(-back_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_3_ID);
		SetAngle(back_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_4_ID);
	//}

}


/*Cyber height setting*/
void CyberUpDown(void)
{

            /*Motor limiter*/
			 if(Up_angle<0)
				 Up_angle=0;
			 else if(Up_angle>1.7)
				 Up_angle=1.7;

     /*Check if Cybergear is ready for control*/
	  // if(Cybercounter==0)
	//	   cntflag=1;



	  // if(cntflag==1)
	   //{
		Cybercounter++;
	   //}


    if(Cybercounter==1)
    {
    	 SetAngle( Up_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_1_ID);
    	 SetAngle(-Up_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_2_ID);

    }else if (Cybercounter==2)
    {
		 SetAngle(-Up_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_3_ID);
		 SetAngle( Up_angle,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_4_ID);

			/*Update values*/
				Up_angle_old=Up_angle;
				front_angle=Up_angle;
				front_angle_old=front_angle;
				back_angle=Up_angle;
				back_angle_old=back_angle;

			    roll_left_old  =back_angle;
			    roll_right_old =back_angle;
			    roll_left      =back_angle;
			  	roll_right     =back_angle;

				cntflag=0;
		 Cybercounter=0;
    }


}

/*Cyber front-angle settings*/
void CyberRollLeft(void)
{

	     if(roll_left>1.7)
		     roll_left=1.7;
		  else if(roll_left<0)
			  roll_left=0.0;


		  SetAngle(-roll_left,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_2_ID);
		  SetAngle(roll_left,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_4_ID);


}


void CyberRollRight(void)
{
          if(roll_right>1.7)
        	  roll_right=1.7;
		  else if(roll_right<0)
			  roll_right=0.0;


	  SetAngle(roll_right,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_1_ID);
	  SetAngle(-roll_right,/*hostID=*/CYBER_HOST_ID, /*motorID=*/CYBER_MOTOR_3_ID);


}



/*Cyber SpeedLimit - softness*/
void SpeedLimit()
{

    /*Check if Cybergear is waiting for control*/
	   if(Cybercounter==0)
		   cntflag=1;



	   if(cntflag==1)
	   {
		Cybercounter++;
	   }


	    if(Cybercounter==1)
	    {

	  	  PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_1_ID);
	  	  PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_2_ID);

	    }else if (Cybercounter==2)
	    {

	      PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_3_ID);
	  	  PositionSpeedLimit(speedLimit, CYBER_HOST_ID, CYBER_MOTOR_4_ID);

				/*Update values*/
	  	     speedLimit_old=speedLimit;
	  	     cntflag=0;
			 Cybercounter=0;
	    }



}


/*Cyber Dummy - just to get real-time data from the motors*/
void CyberDummyRequest()
{

	   Cybercounter++;


	    if(Cybercounter==1)
	    {
	    	CyberFrontAngle();

	    }else if (Cybercounter==2)
	    {
	    	CyberBackAngle();
			Cybercounter=0;
	    }

}
