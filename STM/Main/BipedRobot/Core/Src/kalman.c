#include "main.h"
#include "string.h"
#include "math.h"
#include "mpu6050.h"


extern MPU6050_t MPU6050; //MPU Sensor
#define RAD_TO_DEG 57.295779513082320876798154814105

/*Kalman structures*/
Kalman_t KalmanX = {
    .Q_angle = 0.00001f,    // 0.001f
    .Q_bias = 0.0007f,     //0.003f
    .R_measure = 0.7f,        //0.003f
    .bias=-3.56f        };

Kalman_t KalmanY = {
    .Q_angle = 0.001f,
    .Q_bias = 0.003f,
    .R_measure = 0.03f,
	.bias=1.2f            };

/*******************************************************************************
* Function Name  : Estimate roll and ptich angle from accelerometer
* Description    : MPU6050
* Input          : None
* Output         : None
* Return         : None
*******************************************************************************/
double roll_acc=0,pitch_acc=0;       //Angles calculated directly form ACC vectors X,Y,Z
double dt;
double roll_raw=0, pitch_raw=0; //Corrected values with Kalman filter
double roll_kalman=0, pitch_kalman=0; //Corrected values with Kalman filter

uint32_t timer;

void Acc_ptich_roll_Kalman(MPU6050_t *DataStruct)
{
   //Measure sampling time
   dt = (double)(HAL_GetTick() - timer) / 1000;  //Convert to ms
   timer = HAL_GetTick();


   /*Roll Axis*/
   double roll_sqrt = sqrt( DataStruct->Accel_X_RAW * DataStruct->Accel_X_RAW + DataStruct->Accel_Z_RAW * DataStruct->Accel_Z_RAW);
   if (roll_sqrt != 0.0)
   {
	   roll_acc = atan(DataStruct->Accel_Y_RAW / roll_sqrt) * RAD_TO_DEG;
   }
   else
   {
	   roll_acc = 0.0;
   }



   /*Pitch Axis*/
   double pitch_acc = atan2(-DataStruct->Accel_X_RAW, DataStruct->Accel_Z_RAW) * RAD_TO_DEG;
   if ((pitch_acc < -90 && DataStruct->KalmanAngleY > 90) || (pitch_acc > 90 && DataStruct->KalmanAngleY < -90))
   {
       KalmanY.angle = pitch_acc;
       DataStruct->KalmanAngleY = pitch_acc;
   }
   else
   {

	   double velocity_gyro_Y = DataStruct->Gy;
	   pitch_kalman = Kalman_filter(&KalmanY, pitch_acc, velocity_gyro_Y); //The Kalman filter
	   pitch_raw     = pitch_raw + dt * velocity_gyro_Y;   //No filter - plain integration of the velocity_Y

   }



   //Change direction of the rotation
   if (fabs(DataStruct->KalmanAngleY) > 90)
       DataStruct->Gx = -DataStruct->Gx;


   double velocity_gyro_X = DataStruct->Gx;
   roll_kalman=Kalman_filter(&KalmanX, roll_acc, velocity_gyro_X );
   roll_raw   = roll_raw + dt * velocity_gyro_X;  //No filter - plain integration of the velocity_X
}


/****************************************************************************************************************************
 ********************  Kalman predict ***************************************************************************************

  Process model
  x(k+1) =A x(k) + Bv(k) + w(k)   //w(k)-process nois -> covariance matrix Q
    y(k)= C x(k) + n(k)          //n(k) -sensor  nois -> covariance matrix R

   ___________________________________________________________________________________________

            x(k+1)          =       A         *     x(k)             +     B    *    v(k)
   [angle(k+1)  bias(k+1)]' = [ 1 -Ts, 0   1] * [angle(k)  bias(k)]' + [Ts ,0]  * [Gyro,0]

  ______________________________________________________________________________________________

  //Angle estimate from gyro sensor (integral)
  angle(k) = angle(k-1) + Ts * (gyro_measure - bias);  //Integral; Bias-nonlinearitiy of sensor; Ts- sampling time

*****************************************************************************************************************************/

double Kalman_filter(Kalman_t *Kalman, double angle_acc, double velocity_gyro)
{

    //State space model dx=Ax+Bu (integration of the velocity)
    Kalman->angle = Kalman->angle + dt * (velocity_gyro - Kalman->bias);


    //Prediciton covariance
    Kalman->P[0][0] += dt * (dt * Kalman->P[1][1] - Kalman->P[0][1] - Kalman->P[1][0] + Kalman->Q_angle);
    Kalman->P[0][1] -= dt * Kalman->P[1][1];
    Kalman->P[1][0] -= dt * Kalman->P[1][1];
    Kalman->P[1][1] += dt * Kalman->Q_bias ;

    //Innovation (or pre-fit residual)
    double S = Kalman->P[0][0] + Kalman->R_measure;

    //Optimal Kalamn Gain
    double K[2];
    K[0] = Kalman->P[0][0] / S;
    K[1] = Kalman->P[1][0] / S;


    //Sate Error variable e = xm - xe
    double error = angle_acc - Kalman->angle;

    //State update/correction:  xe = xe + K*error
    Kalman->angle = Kalman->angle + K[0] * error;
    Kalman->bias =  Kalman->bias  + K[1] * error;


    Kalman->velocity = velocity_gyro - Kalman->bias;

    double P00_temp = Kalman->P[0][0];
    double P01_temp = Kalman->P[0][1];

    //Covrainace matrix update
    Kalman->P[0][0] -= K[0] * P00_temp;
    Kalman->P[0][1] -= K[0] * P01_temp;
    Kalman->P[1][0] -= K[1] * P00_temp;
    Kalman->P[1][1] -= K[1] * P01_temp;


    //Output function
    return Kalman->angle;
};
