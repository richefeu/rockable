// Copyright (C) Rockable <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// Rockable can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#include "DrivingSystem.hpp"
#include "Rockable.hpp"

#include "fileTool.hpp"
#include "message.hpp"

DrivingSystem::DrivingSystem() : ServoFunction(nullptr) {}

void DrivingSystem::read(bool warn) {
  static std::map<std::string, int> typeMap{
      {"_x_Vel_", _x_Vel_},       {"_y_Vel_", _y_Vel_},       {"_z_Vel_", _z_Vel_},       {"_xrot_Vel_", _xrot_Vel_},
      {"_yrot_Vel_", _yrot_Vel_}, {"_zrot_Vel_", _zrot_Vel_}, {"_x_For_", _x_For_},       {"_y_For_", _y_For_},
      {"_z_For_", _z_For_},       {"_xrot_Mom_", _xrot_Mom_}, {"_yrot_Mom_", _yrot_Mom_}, {"_zrot_Mom_", _zrot_Mom_},
  };

  if (!fileTool::fileExists("drivingSystem.txt")) {
    if (warn) {
      std::cout << msg::info() << "No file 'drivingSystem.txt' has been found." << std::endl;
      std::cout << " >      The first 'nDriven' bodies are frozen (velocities = 0)." << msg::normal() << std::endl
                << std::endl;
    }
    return;
  }

  std::ifstream is("drivingSystem.txt");

  std::string token;
  is >> token;

  controls.clear();
  while (is) {
    if (token[0] == '/' || token[0] == '#' || token[0] == '!')
      getline(is, token);
    else if (token == "Control") {
      std::string typeStr;
      Control C;
      is >> typeStr >> C.i >> C.value;
      C.type = typeMap[typeStr];
      controls.push_back(C);
    } else if (token == "Servo") {
      std::string servoName;
      is >> servoName;
      if (servoName.substr(0, 6) == "tritri") {
        size_t idXmin, idXmax;
        size_t idYmin, idYmax;
        size_t idZmin, idZmax;
        is >> idXmin >> idXmax >> idYmin >> idYmax >> idZmin >> idZmax;
        Control C;
        C.value = 0.0;  // will be set by the servo
        size_t icontr = controls.size();

        C.i = idXmin;
        C.type = _x_Vel_;
        controls.push_back(C);  // icontr
        C.i = idXmax;
        C.type = _x_For_;
        controls.push_back(C);  // icontr + 1

        C.i = idYmin;
        C.type = _y_Vel_;
        controls.push_back(C);  // icontr + 2
        C.i = idYmax;
        C.type = _y_For_;
        controls.push_back(C);  // icontr + 3

        C.i = idZmin;
        C.type = _z_Vel_;
        controls.push_back(C);  // icontr + 4
        C.i = idZmax;
        C.type = _z_For_;
        controls.push_back(C);  // icontr + 5

        // rotation blocked on the 8 walls
#define _FIXED_ROT(W)    \
  C.i = W;               \
  C.type = _xrot_Vel_;   \
  controls.push_back(C); \
  C.type = _yrot_Vel_;   \
  controls.push_back(C); \
  C.type = _zrot_Vel_;   \
  controls.push_back(C)

        _FIXED_ROT(idXmin);
        _FIXED_ROT(idXmax);
        _FIXED_ROT(idYmin);
        _FIXED_ROT(idYmax);
        _FIXED_ROT(idZmin);
        _FIXED_ROT(idZmax);
#undef _FIXED_ROT

        if (servoName == "tritriIsostaticCompression") {
          // This is an isostatic compression.
          // the pressure is applied on the walls with the maximum positions
          // along the x-, y- and z-axis.
          // The walls with the minimum positions are fixed (velocity imposed at zero)
          double pressure;
          is >> pressure;
          ServoFunction = [icontr, idXmin, idXmax, idYmin, idYmax, idZmin, idZmax, pressure](Rockable& box) -> void {
            const double Xmin = box.Particles[idXmin].pos.x + box.Particles[idXmin].MinskowskiRadius();
            const double Xmax = box.Particles[idXmax].pos.x - box.Particles[idXmax].MinskowskiRadius();
            const double Ymin = box.Particles[idYmin].pos.y + box.Particles[idYmin].MinskowskiRadius();
            const double Ymax = box.Particles[idYmax].pos.y - box.Particles[idYmax].MinskowskiRadius();
            const double Zmin = box.Particles[idZmin].pos.z + box.Particles[idZmin].MinskowskiRadius();
            const double Zmax = box.Particles[idZmax].pos.z - box.Particles[idZmax].MinskowskiRadius();
            const double dX = fabs(Xmax - Xmin);
            const double dY = fabs(Ymax - Ymin);
            const double dZ = fabs(Zmax - Zmin);
            const double surfX = dY * dZ;
            const double surfY = dX * dZ;
            const double surfZ = dY * dX;

            box.System.controls[icontr + 1].value = -pressure * surfX;
            box.System.controls[icontr + 3].value = -pressure * surfY;
            box.System.controls[icontr + 5].value = -pressure * surfZ;
          };
        } else if (servoName == "tritriBiaxialCompression") {
          // This is a biaxial compression along the y-axis
          // the pressure is applied laterally along the x- and z-axis
          // The walls with the minimum positions are fixed (velocity imposed at
          // zero)
          double pressure, velocity;
          is >> pressure >> velocity;
          controls[icontr + 3].type = _y_Vel_;  // because it was previouly set to _y_For_
          controls[icontr + 3].value = -velocity;
          ServoFunction = [icontr, idXmin, idXmax, idYmin, idYmax, idZmin, idZmax, pressure,
                           velocity](Rockable& box) -> void {
            const double Xmin = box.Particles[idXmin].pos.x + box.Particles[idXmin].MinskowskiRadius();
            const double Xmax = box.Particles[idXmax].pos.x - box.Particles[idXmax].MinskowskiRadius();
            const double Ymin = box.Particles[idYmin].pos.y + box.Particles[idYmin].MinskowskiRadius();
            const double Ymax = box.Particles[idYmax].pos.y - box.Particles[idYmax].MinskowskiRadius();
            const double Zmin = box.Particles[idZmin].pos.z + box.Particles[idZmin].MinskowskiRadius();
            const double Zmax = box.Particles[idZmax].pos.z - box.Particles[idZmax].MinskowskiRadius();
            const double dX = fabs(Xmax - Xmin);
            const double dY = fabs(Ymax - Ymin);
            const double dZ = fabs(Zmax - Zmin);
            const double surfX = dY * dZ;
            const double surfZ = dY * dX;

            box.System.controls[icontr + 1].value = -pressure * surfX;
            box.System.controls[icontr + 5].value = -pressure * surfZ;
          };
        } else if (servoName == "tritriCustom") {
          int xminType, xmaxType;
          int yminType, ymaxType;
          int zminType, zmaxType;

          double xminValue, xmaxValue;
          double yminValue, ymaxValue;
          double zminValue, zmaxValue;

          is >> xminType >> xminValue;
          is >> xmaxType >> xmaxValue;
          is >> yminType >> yminValue;
          is >> ymaxType >> ymaxValue;
          is >> zminType >> zminValue;
          is >> zmaxType >> zmaxValue;

          controls[icontr].type = (xminType == 0) ? _x_Vel_ : _x_For_;
          controls[icontr].value = xminValue;
          controls[icontr + 1].type = (xmaxType == 0) ? _x_Vel_ : _x_For_;
          controls[icontr + 1].value = -xmaxValue;
          controls[icontr + 2].type = (yminType == 0) ? _y_Vel_ : _y_For_;
          controls[icontr + 2].value = yminValue;
          controls[icontr + 3].type = (ymaxType == 0) ? _y_Vel_ : _y_For_;
          controls[icontr + 3].value = -ymaxValue;
          controls[icontr + 4].type = (zminType == 0) ? _z_Vel_ : _z_For_;
          controls[icontr + 4].value = zminValue;
          controls[icontr + 5].type = (zmaxType == 0) ? _z_Vel_ : _z_For_;
          controls[icontr + 5].value = -zmaxValue;

          ServoFunction = [icontr, idXmin, idXmax, idYmin, idYmax, idZmin, idZmax, xminValue, xmaxValue, yminValue,
                           ymaxValue, zminValue, zmaxValue](Rockable& box) -> void {
            const double Xmin = box.Particles[idXmin].pos.x + box.Particles[idXmin].MinskowskiRadius();
            const double Xmax = box.Particles[idXmax].pos.x - box.Particles[idXmax].MinskowskiRadius();
            const double Ymin = box.Particles[idYmin].pos.y + box.Particles[idYmin].MinskowskiRadius();
            const double Ymax = box.Particles[idYmax].pos.y - box.Particles[idYmax].MinskowskiRadius();
            const double Zmin = box.Particles[idZmin].pos.z + box.Particles[idZmin].MinskowskiRadius();
            const double Zmax = box.Particles[idZmax].pos.z - box.Particles[idZmax].MinskowskiRadius();
            const double dX = fabs(Xmax - Xmin);
            const double dY = fabs(Ymax - Ymin);
            const double dZ = fabs(Zmax - Zmin);
            const double surfX = dY * dZ;
            const double surfY = dX * dZ;
            const double surfZ = dY * dX;

            if (box.System.controls[icontr].type == _x_For_) box.System.controls[icontr].value = xminValue * surfX;
            if (box.System.controls[icontr + 1].type == _x_For_)
              box.System.controls[icontr + 1].value = -xmaxValue * surfX;

            if (box.System.controls[icontr + 2].type == _y_For_)
              box.System.controls[icontr + 2].value = yminValue * surfY;
            if (box.System.controls[icontr + 3].type == _y_For_)
              box.System.controls[icontr + 3].value = -ymaxValue * surfY;

            if (box.System.controls[icontr + 4].type == _z_For_)
              box.System.controls[icontr + 4].value = zminValue * surfZ;
            if (box.System.controls[icontr + 5].type == _z_For_)
              box.System.controls[icontr + 5].value = -zmaxValue * surfZ;
          };

        } else if (servoName == "tritriLodeAngle") {
          // ADD COMMENTS HERE...
          // **** NOT YET TESTED !!!!!!
          double pressure;
          double LodeAngle;
          double sigRate;
          is >> pressure >> LodeAngle >> sigRate;
          if (LodeAngle < 0.0 || LodeAngle > 60.0) {
            std::cerr << "@DrivingSystem::read, LodeAngle must be set in the range [0° 60°]" << std::endl;
          }
          ServoFunction = [icontr, idXmin, idXmax, idYmin, idYmax, idZmin, idZmax, pressure, LodeAngle, sigRate]
          (Rockable& box) -> void {
            const double Xmin = box.Particles[idXmin].pos.x + box.Particles[idXmin].MinskowskiRadius();
            const double Xmax = box.Particles[idXmax].pos.x - box.Particles[idXmax].MinskowskiRadius();
            const double Ymin = box.Particles[idYmin].pos.y + box.Particles[idYmin].MinskowskiRadius();
            const double Ymax = box.Particles[idYmax].pos.y - box.Particles[idYmax].MinskowskiRadius();
            const double Zmin = box.Particles[idZmin].pos.z + box.Particles[idZmin].MinskowskiRadius();
            const double Zmax = box.Particles[idZmax].pos.z - box.Particles[idZmax].MinskowskiRadius();
            const double dX = fabs(Xmax - Xmin);
            const double dY = fabs(Ymax - Ymin);
            const double dZ = fabs(Zmax - Zmin);
            const double surfX = dY * dZ;
            const double surfY = dX * dZ;
            const double surfZ = dY * dX;
            
            double dSig = sigRate * box.t; // WARNING: initial time MUST be ZERO; dSig(t=0) = 0
            double LodeRadians = LodeAngle * M_PI / 180.0;
            double a = (3.0 * tan(LodeRadians) - sqrt(3.0)) / (2.0 * sqrt(3.0));
            double b = -(1.0 + a);

            box.System.controls[icontr + 1].value = -(pressure + dSig) * surfX;
            box.System.controls[icontr + 3].value = -(pressure + a * dSig) * surfY;
            box.System.controls[icontr + 5].value = -(pressure + b * dSig) * surfZ;
          };
        } else {
          std::cerr << "servoName '" << servoName << "' is unknown" << std::endl;
        }

      } else if (servoName == "shaker") {
        vec3r dir;       // direction of oscillation
        double A, freq;  // amplitude and frequency
        Control C;
        is >> C.i >> dir >> A >> freq;
        dir.normalize();  // in case the user enter a non-normalized one
        double omega = 2.0 * M_PI * freq;
        C.value = 0.0;  // will be set by the servo controller
        size_t icontr = controls.size();
        C.type = _x_Vel_;
        controls.push_back(C);
        C.type = _y_Vel_;
        controls.push_back(C);
        C.type = _z_Vel_;
        controls.push_back(C);
        ServoFunction = [icontr, dir, A, omega](Rockable& box) -> void {
          vec3r v = A * omega * cos(omega * box.t) * dir;
          box.System.controls[icontr].value = v.x;
          box.System.controls[icontr + 1].value = v.y;
          box.System.controls[icontr + 2].value = v.z;
        };
      } else if (servoName == "ramp") {
        std::string typeStr;
        double valueBegin, valueEnd, tBegin, tEnd;
        Control C;
        is >> typeStr >> C.i >> valueBegin >> valueEnd >> tBegin >> tEnd;
        C.type = typeMap[typeStr];
        C.value = valueBegin;
        size_t icontr = controls.size();
        controls.push_back(C);
        ServoFunction = [icontr, valueBegin, valueEnd, tBegin, tEnd](Rockable& box) -> void {
          static double slope = (valueEnd - valueBegin) / (tEnd - tBegin);
          if (box.t <= tBegin)
            box.System.controls[icontr].value = valueBegin;
          else if (box.t >= tEnd)
            box.System.controls[icontr].value = valueEnd;
          else
            box.System.controls[icontr].value = valueBegin + slope * (box.t - tBegin);
        };
      }
    }
    is >> token;
  }
}