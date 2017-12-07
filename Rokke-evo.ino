#include <Servo.h>

int amountServos = 6;
Servo servos[6]; // six servos in array
int curPos[6] = {90, 90, 90}; // stores servo positions

int geneLength = 10; // needs to be an even number
int gene[10];
int servoCenter = 90;
int maxServoRange = 30; // only move maxServoRange from center
int stepSize = 10;

#define PI2 6.283185 // 2*PI saves calculation later

int delayTime = 10; //Delay between each sample when playing back movement (for Tower Pro SG90: At 4.8V, the speed of the servo is 0.12 sec/60° = 2 ms/deg.)
int P = 400; //Maximum period of the periodic movement function
int playBackTime = 2 * P; //The movement pattern will be played back two times when it is being evalutated
#define LENGTH 500 // Length of the wave lookup table - needs to be equal to P/delayTime
byte wave[6][LENGTH]; // Storage for waveform to play back from
float f; //Temporary storage for calculating function values
int offsets[6];

int sensorPin = A0; // select the input pin for LDR (or strain sensor)
int sensorValue = 0; // variable to store the value coming from the sensor
int sensorSum = 0; //Sum of sensor values for calculating fitness

void setup() {
  Serial.begin(9600); //Sets serial port for communication
  for (int i = 0; i < amountServos; i++) {
    servos[i].attach(i + 2);  // Attach the servo on pin i to the servo object
  }
  setServosToPosition(servoCenter);
  delay(100);
}

void loop() {
  int n = 0;
  bool eval = false; 
  while (Serial.available() > 0) {
    delay(20); // ensure entire serial signal is received, not half a signal
    if (n < geneLength) {
      gene[n] = Serial.read();
      Serial.print(gene[n]);
      Serial.print(",");
      eval = true; 
    }
    else {
      Serial.read();
    }
    n++;
  }
  if (eval == true) {
    delay(100);
    float fitness = evaluateIndividual(gene);
    Serial.println(fitness);
    delay(100);
    setServosToPosition(servoCenter);
    delay(100);
    eval = false; 
  }
}

void setServosToPosition(int pos) {
  // smoothly set servo to desired position
  for (int n = 0; n < 10; n++) {
    for (int i = 0; i < amountServos; i++) {
      int dpos = pos - curPos[i];
      int tpos =  pos - ((9 - n) * dpos * 0.11);
      servos[i].write(tpos);
      if (n == 9) {
        curPos[i] = tpos;
      }
    }
    delay(delayTime);
  }
}

float evaluateIndividual(int gene[]) {
  // *1* Wait for and receive output (genotype) from the NN over serial connection
  //  (the format is a 10x2 float array: [wavenumber(1-10)][amplitude[0;maxAngleChange], phase[0;2*pi]])
  //  [copy code here from Franks program to receive values]
  // *2* Calculate array with movement data based on output from NN
  //     (done by a superposition of ten sine waves, i.e. the first ten terms of the Fourier series)
  for (int j = 0; j < amountServos; j++) {
    for (int i = 0; i < LENGTH; i++) { // Step across wave table
      float sum = 0; //Temporary storage for calculating sum of function values
      for (int n = 0; n < geneLength; n += 2) {
        //Standard form of the Fourier series is: sN(x)= A0/2 +Sum(n=1toN)( An*sin(2*pi*n*x/P+phi[n]) );
        //   [N is fixed to 10]
        //   [phi(n) is confined to - [0;2*pi]]
        //   [x is equal to i]
        //   [A0 is set to 90. (Should it be able to vary? Now the fin movements always oscillate around 90 degrees)]
        //   [An is set to vary in the interval between 0 and maxAngleChange]
        // f = gene[n] * sin(((PI2 * (n + 1) * (i + 1)) / P) + gene[n + 1]); //+1 added to n and i to avoid zero for first value
        float AMP = ((float)gene[n] / 256.0) * 30.0;
        float PHASE = 1;//gene[n + 1];
        float FREQ = 1.0/(float)gene[n+1];//gene[n + 2];
        //if (j >=3) {
        f = AMP * sin((PI2 / 2) * ((float) (-i * FREQ)) + PHASE + j); //+1 added to n and i to avoid zero for first value
        //}
        //else {
        //  f = AMP * sin((PI2 / 2) * ((float) (-i * FREQ)));// + PHASE -j); //+1 added to n and i to avoid zero for first value          
        //}
        //Serial.print("AMP:");
        //Serial.print(AMP);
        //Serial.print(",PHASE:");
        //Serial.print(PHASE);
        //Serial.print(",FREQ:");
        //Serial.print(FREQ);
        //Serial.print(",f:");
        //Serial.print(f);
        sum = sum + f;
      }
      sum = 90.0 + (sum / (geneLength / 2.0)) + 0.5; //0.5 added to round off correctly when converting to byte below
      wave[j][i] = (byte)sum;
      //Serial.print("wave ");
      //Serial.print(j);
      //Serial.print(":");
      //Serial.print(wave[j][i]);
    }
  }

  //*3* Play back movement on robot from wavetable and sum up the sensor readings
  for (int i = 0; i < LENGTH; i++) { // Step across wave table
    for (int n = 0; n < amountServos; n++) {
      int servoPos = wave[n][i];
      Serial.print(servoPos);
      Serial.print(",");
      if (servoPos < servoCenter - maxServoRange) {
        servoPos = servoCenter - maxServoRange;
      }
      else if (servoPos < servoCenter - maxServoRange) {
        servoPos = servoCenter - maxServoRange;
      }
      servos[n].write(servoPos); // Move servos to the right position for the current time
      curPos[n] = servoPos;
      sensorValue = analogRead(sensorPin); // Read the value from the sensor (comes in the range of 0 to 1023)
      // Serial.println(sensorValue); //Print the values coming from the sensor on the screen
      sensorSum = sensorSum + sensorValue; //Add sensor reading to the sum (summed here because saving it clogged up memory)
      delay(delayTime); //Wait before playing the next sample [Is this needed? How long does it take to execute the code in the loop? Maybe good to have to protect the servos?]
    }
  }

  //*4* Send the summation of sensor readings to NN (this value is used to calculate the fitness of the individual (the phenotype))
  // [copy code from Franks program]
  return sensorSum = 0;
}

