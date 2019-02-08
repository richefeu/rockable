###This folder contains header-only c++ tools used by most of the codes in the folder mbox.

**AABB.hpp**: an Axis Aligned Bounding Box

**ColorTable.hpp**: 

**ExecChrono.hpp**:

```
 ExecChrono MM; 
 // the chrono will start automatically (but we can use 'start' if we want)
 // Doing something...
 MM.stop();

```

**Mth.hpp**: some math constants

**OBB.hpp**: an Oriented Bounding Box

**any_type.hpp**:

**DataTable.hpp**:

```
DataTable dt;

dt.set("mu", 0, 0, 0.7);
dt.set("kn", 0, 1, 1000);
dt.set("kn", 1, 2, 200.8);

dt.write(std::cout);
```

**factory.hpp**: 

**fastSort3.hpp**: fast sort of 3 values (template)

```
 int a = 3, b = 1, c = 2;
 fastSort3<int>(a, b, c);
 // result: a = 1, b = 2, c = 3
```

**fileTool.hpp**: file tools

**geoTool.hpp**: some geometric tool 

**grid2.hpp**: nxn grid

**histo.hpp**: compute probability density function or probability distribution

```
std::vector<double> v;
for (size_t i = 0 ; i < 1000 ; i++) v.push_back(i/1000.);

histo H = pdfMaxPerBin(v, 200);
for (size_t i = 0 ; i < H.P.size() ; i++) std::cout << H.P[i].X << " " << H.P[i].P << " " << H.P[i].W << std::endl;
```

**linreg.hpp**: linear regression

**mat4.hpp**: matrix 2x2 (double)

**mat4sym.hpp**: symmetric matrix 2x2

**mat9.hpp**: matrix 3x3 (template)

**mat9sym.hpp**: symmetric matrix 3x3 (template)

**message.hpp**: some tools for printing messages with colors in the terminal

**nextToken.hpp**:

```
std::string token;
nextToken tokenizer(stream);
tokenizer.getNext(token);
```

**psio.hpp**:

**quat.hpp**:

**size_generator.hpp**:

**slicedRange.hpp**:

**stackTracer.hpp**:

**triangulatePolygon.hpp**:

**typedefs.hpp** (deprecated)

**util.hpp**:

**vec2.hpp**:

**vec3.hpp**:

