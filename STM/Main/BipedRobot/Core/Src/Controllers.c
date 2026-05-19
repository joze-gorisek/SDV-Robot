/* STABILIZACIJA CONTROLLER    */


#include "Controllers.h"
#include "DDSM115.h"
#include "string.h"
#include "mpu6050.h"
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

extern double roll_kalman, pitch_kalman;
extern double roll_acc, pitch_acc;
extern Kalman_t KalmanX;
extern Kalman_t KalmanY;


/*Drive Wheels*/
extern float DDSangle01;       //Current angle    [0-32767] -> [0 360]   address= 0x01 over RS485
extern float DDSrpm01;         //Current velocity [0-     ] -> [-330RPM  330RPM]
extern float DDScurrent01;     //Current Torque   [-32767   32767] -> [-8A   8A]
extern float DDSvelocityRadial01;    //Velocity
extern float distnacex01;
extern float DDSvelocity01;

extern float DDSangle30;       //Current angle    [0-32767] -> [0 360]   address= 0x30 over RS485
extern float DDSrpm30;         //Current velocity [0-     ] -> [-330RPM  330RPM]
extern float DDScurrent30;     //Current Torque   [-32767   32767] -> [-8A   8A]
extern float DDSvelocityRadial30;    //Velocity
extern float distnacex30;
extern float DDSvelocity30;

int ref_flag=0;
double iq = 0,iq01 = 0,iq30 = 0;
double iqmot = 0;
double K[2]   ={13.0053,  2.0984}; //LQR state space
double Kack[2]={7.2788,    0.1049}; //LQR state space 7.2307,    0.0092
double offset=-1.5;
double sens=1.0;
double velocity_k   = 0;
double velocity_k_1 = 0;
double pole_filter = -160;
double Ts = 0.020;
double velocity_kalman=0;
double velocity_kalman_1=0.0;
double angle_robot=0.0f;
double integral_robust=0.0f,integral_robust01=0.0,integral_robust30=0.0;
double integral_x01=0.0,integral_x30=0.0f;
double velocity_robot=0;
double Coffset=0.1;     //0.123
double DeadZone=0.7;  //0.1
double angle_old=0;
double SSref01=0.0f,SSref30=0.0f;  //Reference value for state-space controller
double SSref_old01=0.0f,SSref_old30=0.0f;
double SSrefd01=0.0f,SSrefd30=0.0f;  //Derivative value of the reference signal
double  Xref=0.0;
double  Xref01=0.0; //Right motor
double Xref30=0.0; //Left  motor
double  Xerror01=0.0,Xerror30=0.0;
double

Xerror=0.0;

double k1=2.1f,k2=0.44f, k3=0.341f, k4=0.00; //Controller gain [k1 k2 k3] 1.5 0.77 0.45  Coffset=0.11
float k1l=0.055f,k2l=0.0000f, k3l=0.00f; //Controller gain [k1 k2 k3] 1.5 0.77 0.45  Coffset=0.11
double k1x=0.082,k2x=0.012,k3x=0.0;

extern double dt;
extern int ControllerON;
extern float Up_angle;
extern float roll_right;
extern float roll_left;
float lateral_out=0.0,error_lateral=0.0,integral_lateral=0.0;
float KalmanY_angle_old=0.0;
float KalmanY_angle_offset=2.89;

void LateralController(uint8_t MOTOR_ID, uint8_t ControllerStatus)
{




}

void LQR_controllerLR(uint8_t MOTOR_ID, uint8_t ControllerStatus) {



}




/*OLD controller*/

double integral_x=0.0,SSref=0.0,SSrefd=0,SSref_old=0;

void LQR_controller(uint8_t MOTOR_ID, uint8_t ControllerStatus) {

	/* LOW PASS FILTER - VELOCITY */
	velocity_k = exp((double) pole_filter * dt ) *velocity_k_1 + (1-exp((double) pole_filter * dt )) * velocity_kalman_1;

			/* Convolution */
			velocity_kalman = KalmanX.velocity; //Filter input
			velocity_kalman_1=velocity_kalman;
			velocity_k_1 = velocity_k;



}






