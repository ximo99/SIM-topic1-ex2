// ----------------------------------------------- //<>//
// Plantilla para la implementación de integradores 
// numéricos en Simulación. Tema 1.
// 
// Sergio Casas & Miguel Lozano (Febrero 2021)
// --------------------------------------------------


// Definitions:
enum IntegratorType 
{
  NONE,
  EXPLICIT_EULER,    //Tiempo:  3.0099976 Error:  2.463363
  SIMPLECTIC_EULER,  //Tiempo:  3.0099976 Error:  2.2371943
  HEUN,              //Tiempo:  3.0099976 Error:  76.78063
  RK2,               //Tiempo:  3.0099976 Error:  2.3450334
  RK4                //Tiempo:  3.0099976 Error:  2.3450265

}

// Parameters of the numerical integration:
float SIM_STEP = 0.01;   // Simulation time-step (s)
IntegratorType integrator = IntegratorType.RK4;   // ODE integration method

// Display values:
final boolean FULL_SCREEN = false;
final int DRAW_FREQ = 50;   // Draw frequency (Hz or Frame-per-second)
int DISPLAY_SIZE_X = 1000;  // Display width (pixels)
int DISPLAY_SIZE_Y = 600;   // Display height (pixels)

// Draw values:
final int [] BACKGROUND_COLOR = {200, 200, 255};
final int [] REFERENCE_COLOR = {0, 255, 0};
final int [] OBJECTS_COLOR = {255, 0, 0};

final float OBJECTS_SIZE = 1.0;   // Size of the objects (m)
final float PIXELS_PER_METER = 20.0;   // Display length that corresponds with 1 meter (pixels)
final PVector DISPLAY_CENTER = new PVector(0.0, 0.0);   // World position that corresponds with the center of the display (m)

//Time control
float _simTime = 0.0;
float _elapsedTime = 0.0; //Elapsed (real) time  (s)

// Parameters of the problem:
final float M   = 10.0;   // Particle mass (kg)
final float Gc  = 9.801;   // Gravity constant (m/(s*s))
final float K   = 0.2;
final PVector G = new PVector(0.0, -Gc);   // Acceleration due to gravity (m/(s*s))

PVector _s = new PVector();   // Position of the particle (m)
PVector _v = new PVector();   // Velocity of the particle (m/s)
PVector _a = new PVector();   // Accleration of the particle (m/(s*s))

final PVector s0 = new PVector(-20.0, -20.0);   // Particle's start position (m)
PVector _v0 = new PVector(8.0, 20.0, 0.0);      // Velocidad inicial
PVector s_a = new PVector();                    // Posición analítica


//------------------------------------------------------
// Inicialización

void settings()
{
  if (FULL_SCREEN)
  {
    fullScreen();
    DISPLAY_SIZE_X = displayWidth;
    DISPLAY_SIZE_Y = displayHeight;
  } 
  else
    size(DISPLAY_SIZE_X, DISPLAY_SIZE_Y);
}

void setup()
{
  frameRate(DRAW_FREQ);
  
  initSimulation(); 
}

//-------------------------------------------------------
// Simulation

void initSimulation()
{
  _simTime = 0.0;
  _elapsedTime = 0.0;
  
  _v0.set(8.0, 20.0);
  _s = s0.copy();
  
  _v.set(_v0.x, _v0.y);
  _a.set(0.0, 0.0);
}

void updateSimulation()
{
  switch (integrator)
  {
    case EXPLICIT_EULER:
      updateSimulationExplicitEuler();
      break;
  
    case SIMPLECTIC_EULER:
      updateSimulationSimplecticEuler();
      break;
  
    case HEUN:
      updateSimulationHeun();
      break;
  
    case RK2:
      updateSimulationRK2();
      break;
  
    case RK4:
      updateSimulationRK4();
      break;
  }
  
  _simTime += SIM_STEP;
  
  // Solucion analítica la posicion: doble integracion
  s_a.x = s0.x + (M * _v0.x / K) * (1-exp((-K / M) * _simTime));
  s_a.y = s0.y + (M/K) * ((M*Gc/K) + _v0.y)*(1 - exp((-K / M) * _simTime)) - (M * Gc * _simTime) / K;
  
  println("Tiempo: ", str(_simTime), "Error: ", PVector.sub(s_a, _s).mag());
  
  if (_simTime > 10.0)
    exit();
}

void updateSimulationExplicitEuler()
{
  // s(t+h) = s(t) + h*v(t)
  // v(t+h) = v(t) + h*a(s(t),v(t))
  
  _a = calculateAcceleration(_s, _v);
  
  _s.add(PVector.mult(_v, SIM_STEP));
  _v.add(PVector.mult(_a, SIM_STEP));
  
}

void updateSimulationSimplecticEuler()
{
  // v(t+h) = v(t) + h*a(s(t),v(t))
  // s(t+h) = s(t) + h*v(t+h)
  
  //Calcular la derivada en el principio del intervalo
  _a = calculateAcceleration(_s, _v);
  
  // Calcular la velocidad siguiente a partir de la derivada en el principio del intervalo
  _v.add(PVector.mult(_a, SIM_STEP));
  
  // Calcular la posición siguiente a partir de la velocidad en el principio del intervalo
  _s.add(PVector.mult(_v, SIM_STEP));
}

void updateSimulationHeun()
{
  //Integración numérica de la velocidad
    //Calcular acceleracion _a
    _a = calculateAcceleration(_s, _v);
    
    //Paso de Euler, actualizo s2, v2
    PVector _s2 = new PVector();
    _s2 = _s;
    _s2.add(PVector.mult(_v, SIM_STEP));
    PVector _v2 = new PVector();
    _v2 = _v;
    _v2.add(PVector.mult(_a, SIM_STEP));
    
    //v_promedio = (_v + v2)/2
    PVector v_prom = PVector.mult(PVector.add(_v, _v2), 0.5);
    
    //actualizar _s con la v promedio
    _s.add(PVector.mult(v_prom, SIM_STEP));
    
  //Integración de la acceleración
    //Calcular la acceleración al final del intervalo
    //a2 = a(s2, v2)
    PVector _a2 = new PVector();
    _a2 = calculateAcceleration(_s2, _v2);
    
    //a_promedio = (_a + a2)/2
    PVector a_prom = PVector.mult(PVector.add(_a, _a2), 0.5);
    
    //Actualizar la velocidad (_v) con la aceleración promedia
    _v.add(PVector.mult(a_prom, SIM_STEP));
}

void updateSimulationRK2()
{
  //Calcular acceleracion _a
  _a = calculateAcceleration(_s, _v);

  // k1v = a(s(t), v(t)) * h
  PVector k1s = PVector.mult(_v, SIM_STEP);
  
  // k1s = v(t) * h
  PVector k1v = PVector.mult(_a, SIM_STEP);
  
  PVector s2 = PVector.add(_s, PVector.mult(k1s, 0.5));
  PVector v2 = PVector.add(_v, PVector.mult(k1v, 0.5));
  PVector a2 = calculateAcceleration(s2, v2);

  // k2v = a(s(t) + k1s / 2, v(t) + k1v / 2) * h
  PVector k2v = PVector.mult(a2, SIM_STEP);
  
  // k2s = (v(t)+k1v/2)*h
  PVector k2s = PVector.mult(PVector.add(_v, PVector.mult(k1v, 0.5)), SIM_STEP);
  
  _v.add(k2v);
  _s.add(k2s);

}

void updateSimulationRK4()
{
  //Calcular acceleracion _a
  _a = calculateAcceleration(_s, _v);

  // k1v = a(s(t), v(t)) * h
  PVector k1v = PVector.mult(_a, SIM_STEP);
  
  // k1s = v(t) * h
  PVector k1s = PVector.mult(_v, SIM_STEP);
  
  PVector s2 = PVector.add(_s, PVector.mult(k1s, 0.5));
  PVector v2 = PVector.add(_v, PVector.mult(k1v, 0.5));
  PVector a2 = calculateAcceleration(s2, v2);

  // k2v = a(s(t) + k1s / 2, v(t) + k1v / 2) * h
  PVector k2v = PVector.mult(a2, SIM_STEP);
  
  // k2s = (v(t)+k1v/2)*h
  PVector k2s = PVector.mult(PVector.add(_v, PVector.mult(k1v, 0.5)), SIM_STEP);
  
  PVector s3 = PVector.add(_s, PVector.mult(k2s, 0.5));
  PVector v3 = PVector.add(_v, PVector.mult(k2v, 0.5));
  PVector a3 = calculateAcceleration(s3, v3);

  // k3v = a(s(t)+k2s/2, v(t)+k2v/2)*h
  PVector k3v = PVector.mult(a3, SIM_STEP);
  
  // k3s = (v(t)+k2v/2)*h
  PVector k3s = PVector.mult(PVector.add(_v, PVector.mult(k2v, 0.5)), SIM_STEP);

  PVector s4 = PVector.add(_s, k3s);
  PVector v4 = PVector.add(_v, k3v);
  PVector a4 = calculateAcceleration(s4, v4);

  // k4v = a(s(t)+k3s, v(t)+k3v)*h
  PVector k4v = PVector.mult(a4, SIM_STEP);
  
  // k4s = (v(t)+k3v)*h
  PVector k4s = PVector.mult(PVector.add(_v,k3v), SIM_STEP);
  
  // v(t+h) = v(t) + (1/6)*k1v + (1/3)*k2v + (1/3)*k3v +(1/6)*k4v
  // s(t+h) = s(t) + (1/6)*k1s + (1/3)*k2s + (1/3)*k3s +(1/6)*k4s
  
  _v.add(PVector.mult(k1v, 1/6.0));
  _v.add(PVector.mult(k2v, 1/3.0));
  _v.add(PVector.mult(k3v, 1/3.0));
  _v.add(PVector.mult(k4v, 1/6.0));
  
  _s.add(PVector.mult(k1s, 1/6.0));
  _s.add(PVector.mult(k2s, 1/3.0));
  _s.add(PVector.mult(k3s, 1/3.0));
  _s.add(PVector.mult(k4s, 1/6.0)); 
}


//-----------------------------------------------------------
// Aceleración del problema dy/dx = dv/dt, función a integrar

// Problema - Tiro arabólico

// Ecuacion Diferencial:
// s' = v(t)
// v' = a(s(t), v(t))
// siendo: 
//      a(s(t), v(t)) = [Froz(v(t)) + Fpeso ]/m
//      Froz = -k·v(t)
//      Fpeso = mg; siendo g(0, -9.8) m/s²

PVector calculateAcceleration(PVector s, PVector v)
{
  PVector Froz   = PVector.mult(v,K);
  PVector Fpeso  = PVector.mult(G, M); 

  PVector f = PVector.add(Fpeso, Froz);
  
  PVector a = PVector.div(f, M);
  return a;
}

//-------------------------------------------------------
void draw()
{
  background(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2]);
  
  updateSimulation();

  drawScene();
}

void drawScene()
{
  fill(OBJECTS_COLOR[0], OBJECTS_COLOR[1], OBJECTS_COLOR[2]);
  strokeWeight(1);

  PVector screenPos = new PVector();
  worldToScreen(_s, screenPos);
  circle(screenPos.x, screenPos.y, worldToPixels(OBJECTS_SIZE));
  
}

//-----------------------------------------------------
// World2Screen functions ...

// Converts distances from world length to pixel length
float worldToPixels(float dist)
{
  return dist*PIXELS_PER_METER;
}

// Converts distances from pixel length to world length
float pixelsToWorld(float dist)
{
  return dist/PIXELS_PER_METER;
}

// Converts a point from world coordinates to screen coordinates
void worldToScreen(PVector worldPos, PVector screenPos)
{
  screenPos.x = 0.5*DISPLAY_SIZE_X + (worldPos.x - DISPLAY_CENTER.x)*PIXELS_PER_METER;
  screenPos.y = 0.5*DISPLAY_SIZE_Y - (worldPos.y - DISPLAY_CENTER.y)*PIXELS_PER_METER;
}

// Converts a point from screen coordinates to world coordinates
void screenToWorld(PVector screenPos, PVector worldPos)
{
  worldPos.x = ((screenPos.x - 0.5*DISPLAY_SIZE_X)/PIXELS_PER_METER) + DISPLAY_CENTER.x;
  worldPos.y = ((0.5*DISPLAY_SIZE_Y - screenPos.y)/PIXELS_PER_METER) + DISPLAY_CENTER.y;
}

//-----------------------------------------------------
// Interaction

void keyPressed()
{
  if (key == 'r' || key == 'R')
  {
    initSimulation();
  }
 
  if (key == '+')
  {
    SIM_STEP *= 1.1;
  }
  if (key == '-')
  {
    SIM_STEP /= 1.1;
  }
}
