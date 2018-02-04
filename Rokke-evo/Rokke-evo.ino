#include <Servo.h>

int amountServos = 6;
Servo servos[6]; // six servos in array
int curPos[6] = {90, 90, 90}; // stores servo positions
int counter = 0;
bool endSwitch = false;
bool startSwitch = false;
const int startSwitchPin = 14;
const int endSwitchPin = 16;
const int trigPin = 44;
const int echoPin = 42;
long duration;
int distance;

#define IRsensor A0

int offset[6] = { -8, -2, 7, 0, -1, -2};

#define geneLength 15 // needs to be an even number
int gene[geneLength];
int servoCenter = 90;
int maxServoRange = 40; // only move maxServoRange from center
int stepSize = 10;
int maxIRValue = 245;
int maxUSValue = 5000;
int pauseTime = 10000;

#define PI2 6.283185 // 2*PI saves calculation later

int delayTime = 100; //Delay between each sample when playing back movement (for Tower Pro SG90: At 4.8V, the speed of the servo is 0.12 sec/60° = 2 ms/deg.)
int P = 400; //Maximum period of the periodic movement function
int playBackTime = 2 * P; //The movement pattern will be played back two times when it is being evalutated
#define LENGTH 200 // Length of the wave lookup table - needs to be equal to P/delayTime
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
  pinMode(startSwitchPin, INPUT_PULLUP);
  pinMode(endSwitchPin, INPUT_PULLUP);
  pinMode(trigPin, OUTPUT); // Sets the trigPin as an Output
  pinMode(echoPin, INPUT); // Sets the echoPin as an Input
  // IMPORTANT FOR EVOLUTION::
  while (!Serial) {
    ; // wait for serial port to connect. Needed for native USB port only
  }
  establishContact();
}

void establishContact() {
  while (Serial.available() <= 0) {
    Serial.println("0,0,0");   // send an initial string
    delay(300);
  }
}

float measureUltrasound() { // has a 8 ms delay
  // Clears the trigPin
  digitalWrite(trigPin, LOW);
  delayMicroseconds(2);
  // Sets the trigPin on HIGH state for 10 micro seconds
  digitalWrite(trigPin, HIGH);
  delayMicroseconds(6);
  digitalWrite(trigPin, LOW);
  // Reads the echoPin, returns the sound wave travel time in microseconds
  duration = pulseIn(echoPin, HIGH);
  // Calculating the distance
  // distance = duration * 0.034 / 2;
  // Prints the distance on the Serial Monitor
  // Serial.print("Distance: ");
  // Serial.println(distance);
  return duration;
}

void loop() {
  int n = 0;
  bool eval = false;
  // Serial.println(analogRead(IRsensor));  // value from sensor * (5/1024)
  // delay(500);
  //Serial.println(measureUltrasound());
  // setServosToPosition(servoCenter);

  while (Serial.available()) {
    delay(500); // ensure entire serial signal is received, not half a signal
    if (n < geneLength) {
      gene[n] = Serial.read();
      // Serial.print(gene[n]);
      // Serial.println(",");
      eval = true;
    }
    else {
      Serial.read();
    }
    n++;
  }
  //eval = true;
  if (eval == true) {
    delay(100);
    Serial.flush();
    moveToStart();
    setServosToPosition(servoCenter);
    delay(pauseTime);
    //  float AMP = ((float)gene[n] / 255.0) * maxServoRange;
    //  float PHASE = gene[n + 1] / 255.0 * 4;
    //  float FREQ = (float)gene[n + 2] / 255; //gene[n + 2];
    float fitness = evaluateIndividual(gene);
    Serial.println(fitness);
    delay(100);
    setServosToPosition(servoCenter);
    delay(100);
    eval = false;
  }
  //
}

void moveToStart() {
  float AMP = 40.0;
  float PHASE = 1.0;//gene[n + 1]; // max 4
  float FREQ = 100.0 / 255.0 ; //1.0;//gene[n + 2]; // max 1
  //  Serial.print("STARTING");
  int i = 0;
  while (startSwitch == false) {
    i++;
    // old // timeCount += FREQ;
    //    Serial.print(digitalRead(startSwitchPin));
    if (digitalRead(startSwitchPin) == LOW) {
      startSwitch = true;
    }
    for (int n = 0; n < amountServos; n++) {
      int f = AMP * sin((PI2 / 2) * ((float) (i * FREQ)) + (PHASE * n)); //+1 added to n and i to avoid zero for first value

      int servoPos = 90 + f;
      // old // int servoPos = 90 + (int)(AMP * sin(timeCount + (n * 1.0)));
      //Serial.print("tc = ");
      //Serial.print(timeCount);
      //Serial.print(servoPos);
      //Serial.print(",");
      if (servoPos < servoCenter - maxServoRange) {
        servoPos = servoCenter - maxServoRange;
        //        Serial.print("--");
      }
      else if (servoPos > servoCenter + maxServoRange) {
        servoPos = servoCenter + maxServoRange;
      }
      //      if (n == 4) {
      servos[n].write(servoPos + offset[n]); // Move servos to the right position for the current time
      curPos[n] = servoPos;
      //      }
      // sensorValue = analogRead(sensorPin); // Read the value from the sensor (comes in the range of 0 to 1023)
      // Serial.println(sensorValue); //Print the values coming from the sensor on the screen
      // sensorSum = sensorSum + sensorValue; //Add sensor reading to the sum (summed here because saving it clogged up memory)
    }
    delay(delayTime); //Wait before playing the next sample [Is this needed? How long does it take to execute the code in the loop? Maybe good to have to protect the servos?]
    //delay(8);
  }
  //  Serial.print("DONE");
  startSwitch = false;
}

void setServosToPosition(int pos) {
  // smoothly set servo to desired position
  for (int i = 0; i < amountServos; i++) {
    servos[i].write(90 + offset[i]);
  }
  delay(delayTime);
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
      for (int n = 0; n < geneLength; n += 3) {
        //Standard form of the Fourier series is: sN(x)= A0/2 +Sum(n=1toN)( An*sin(2*pi*n*x/P+phi[n]) );
        //   [N is fixed to 10]
        //   [phi(n) is confined to - [0;2*pi]]
        //   [x is equal to i]
        //   [A0 is set to 90. (Should it be able to vary? Now the fin movements always oscillate around 90 degrees)]
        //   [An is set to vary in the interval between 0 and maxAngleChange]
        // f = gene[n] * sin(((PI2 * (n + 1) * (i + 1)) / P) + gene[n + 1]); //+1 added to n and i to avoid zero for first value
        float AMP = ((float)gene[n] / 255.0) * maxServoRange;
        float PHASE = gene[n + 1] / 255.0 * 4;
        float FREQ = (float)gene[n + 2] / 255; //gene[n + 2];
        //if (j >=3) {
        f = AMP * sin((PI2 / 2) * ((float) (-i * FREQ)) + (PHASE * j)); //+1 added to n and i to avoid zero for first value
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
        //Serial.println(f);
        if (geneLength <= 3) {
          sum = sum + (f / (geneLength / 3));
        }
        else {
          sum = sum + (f / (geneLength / 3) * 2);
        }
      }
      sum = 90.0 + (sum + 0.5); //0.5 added to round off correctly when converting to byte below
      wave[j][i] = (byte)sum;
      //Serial.print("wave ");
      //Serial.print(j);
      //Serial.print(":");
      //Serial.print(wave[j][i]);
    }
  }

  //*3* Play back movement on robot from wavetable and sum up the sensor readings
  int lCount;
  double fitValue = 0;
  //Serial.print("STARTING EVALUATION");
  for (int i = 0; i < LENGTH; i++) { // Step across wave table
    //Serial.print(digitalRead(endSwitchPin));
    if (digitalRead(endSwitchPin) == LOW) {
      return ((double)(LENGTH - i) * (double) maxUSValue) + fitValue;
      break;
    }
    else {
      lCount++;
      // fitValue += analogRead(IRsensor);   // print the distance
      fitValue += measureUltrasound();
      delay(delayTime);
    }
    for (int n = 0; n < amountServos; n++) {
      int servoPos = wave[n][i];
      //      Serial.println(servoPos);
      //      Serial.print(",");
      if (servoPos < servoCenter - maxServoRange) {
        servoPos = servoCenter - maxServoRange;
      }
      else if (servoPos > servoCenter + maxServoRange) {
        servoPos = servoCenter + maxServoRange;
      }
      servos[n].write(servoPos  + offset[n]); // Move servos to the right position for the current time
      curPos[n] = servoPos;
    }
    // sensorValue = analogRead(sensorPin); // Read the value from the sensor (comes in the range of 0 to 1023)
    // Serial.println(sensorValue); //Print the values coming from the sensor on the screen
    // sensorSum = sensorSum + sensorValue; //Add sensor reading to the sum (summed here because saving it clogged up memory)
    // delay(delayTime); //Wait before playing the next sample [Is this needed? How long does it take to execute the code in the loop? Maybe good to have to protect the servos?]
  }
  return fitValue;
  return 0.0;
  //*4* Send the summation of sensor readings to NN (this value is used to calculate the fitness of the individual (the phenotype))
  // [copy code from Franks program]
  return sensorSum = 0;
}

