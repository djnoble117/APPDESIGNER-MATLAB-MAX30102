/*
  MAX30105 Breakout: Output all the raw Red/IR/Green readings
  By: Nathan Seidle @ SparkFun Electronics
  Date: October 2nd, 2016
  https://github.com/sparkfun/MAX30105_Breakout

  Outputs all Red/IR/Green values.

  Hardware Connections (Breakoutboard to Arduino):
  -5V = 5V (3.3V is allowed)
  -GND = GND
  -SDA = A4 (or SDA)
  -SCL = A5 (or SCL)
  -INT = Not connected

  The MAX30105 Breakout can handle 5V or 3.3V I2C logic. We recommend powering the board with 5V
  but it will also run at 3.3V.

  This code is released under the [MIT License](http://opensource.org/licenses/MIT).
*/

#include <Wire.h>
#include "MAX30105.h"

MAX30105 particleSensor;

#define debug Serial //Uncomment this line if you're using an Uno or ESP
//#define debug SerialUSB //Uncomment this line if you're using a SAMD21

void setup()
{
 Serial.begin(115200);
 // Initialize sensor
 if (!particleSensor.begin(Wire, I2C_SPEED_FAST))
//Use default I2C port, 400kHz speed

  // Initialize sensor
  if (particleSensor.begin() == false)
  {
    debug.println("MAX30105 was not found. Please check wiring/power. ");
    while (1);
  }

  //Setup to sense a nice looking saw tooth on the plotter
 byte ledBrightness = 250; //Options: 0=Off to 255=50mA
 byte sampleAverage = 1; //Options: 1, 2, 4, 8, 16, 32
 byte ledMode = 2; //Options: 1, 2, 3
 int sampleRate = 400; //Options: 50, 100, 200, 400, 800, 1000, 1600
 int pulseWidth = 69; //Options: 69, 118, 215, 411
 int adcRange = 16384; //Options: 2048, 4096, 8192, 16384
 particleSensor.setup(ledBrightness, sampleAverage, ledMode, sampleRate, pulseWidth,
adcRange); //Configure sensor with these settings
}

void loop()
{
  debug.print(" R[");
  debug.print(particleSensor.getRed());
  debug.print("] IR[");
  debug.print(particleSensor.getIR());
  debug.print("]");
  debug.println();
}
