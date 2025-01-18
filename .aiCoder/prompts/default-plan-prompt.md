```markdown
# Plan for a NURBS Library in JavaScript

## Introduction
This document outlines the classes and methods required for a complete NURBS (Non-Uniform Rational B-Spline) library in JavaScript, designed to be the foundation for a Boundary Representation (BREP) kernel. This plan is inspired by the structure of libraries like OpenCASCADE, aiming for a robust and feature-rich implementation.

## Core Concepts

### 1. Control Points
- **Description**: Points in 3D space that define the shape of the NURBS curve or surface.
- **Data**: Typically stored as an array of 3D coordinates (x, y, z), optionally with a weight.

### 2. Knots
- **Description**: A sequence of numbers that define the parameter space for the NURBS curve or surface.
- **Data**:  An array of numbers that are monotonically increasing.

### 3. Weights
- **Description**: Values associated with control points, influencing the curve/surface shape.
- **Data**: An array of numbers, one for each control point.

### 4. Degree
- **Description**: The polynomial degree of the B-spline basis functions.
- **Data**: Integer.

## Class Outline

### 1.  `Point3D` Class
   - **Description**: Represents a point in 3D space.
   - **Methods**:
        - `constructor(x, y, z)`: Creates a new point object with given coordinates.
            -   `x` (Number) : X coordinate of the point.
            -   `y` (Number) : Y coordinate of the point.
            -   `z` (Number) : Z coordinate of the point.
        - `getX()`: Returns the x coordinate of the point.
        - `getY()`: Returns the y coordinate of the point.
        - `getZ()`: Returns the z coordinate of the point.
        - `setX(x)`: Sets the x coordinate of the point.
            -   `x` (Number) : New x coordinate of the point.
        - `setY(y)`: Sets the y coordinate of the point.
            -   `y` (Number) : New y coordinate of the point.
        - `setZ(z)`: Sets the z coordinate of the point.
            -   `z` (Number) : New z coordinate of the point.
        - `toArray()`: Returns the point as an array `[x,y,z]`.
        - `distanceTo(otherPoint)`: Calculate the distance to another point
            - `otherPoint` (Point3D) :  The other point to calculate distance to.
            - Returns: `Number` - distance between the two points.
        - `translate(vector)`: Returns a new translated point
            - `vector` (Vector3D) :  Translation vector.
            - Returns: `Point3D` - new translated point
        - `add(vector)`: Returns a new translated point
             - `vector` (Vector3D) :  Translation vector.
             - Returns: `Point3D` - new translated point
        -  `sub(vector)`: Returns a new translated point
             - `vector` (Vector3D) :  Translation vector.
             - Returns: `Point3D` - new translated point

### 2. `Vector3D` Class
   - **Description**: Represents a vector in 3D space.
    - **Methods**:
        - `constructor(x, y, z)`: Creates a new vector object with given coordinates.
            -   `x` (Number) : X coordinate of the vector.
            -   `y` (Number) : Y coordinate of the vector.
            -   `z` (Number) : Z coordinate of the vector.
        - `getX()`: Returns the x coordinate of the vector.
        - `getY()`: Returns the y coordinate of the vector.
        - `getZ()`: Returns the z coordinate of the vector.
        - `setX(x)`: Sets the x coordinate of the vector.
            -   `x` (Number) : New x coordinate of the vector.
        - `setY(y)`: Sets the y coordinate of the vector.
            -   `y` (Number) : New y coordinate of the vector.
        - `setZ(z)`: Sets the z coordinate of the vector.
            -   `z` (Number) : New z coordinate of the vector.
        - `toArray()`: Returns the vector as an array `[x,y,z]`.
        - `magnitude()`: Calculate the magnitude of the vector.
            - Returns: `Number` - magnitude of the vector.
         - `normalize()`: Normalizes the vector to unit length.
            - Returns: `Vector3D` - normalized vector
        - `dot(otherVector)`: Returns the dot product with another vector
             - `otherVector` (Vector3D) :  The other vector to use for dot product.
             - Returns: `Number` - the dot product result
        - `cross(otherVector)`: Returns the cross product with another vector
             - `otherVector` (Vector3D) :  The other vector to use for cross product.
             - Returns: `Vector3D` - the cross product result
        - `add(otherVector)`: Returns the sum of this vector and another vector
             - `otherVector` (Vector3D) :  The other vector to use for addition.
             - Returns: `Vector3D` - the result of the addition.
       - `sub(otherVector)`: Returns the difference of this vector and another vector
            - `otherVector` (Vector3D) :  The other vector to use for subtraction.
            - Returns: `Vector3D` - the result of the subtraction.
       -  `scale(scalar)`: Scales the vector by a scalar value
             - `scalar` (Number) :  The value to scale by.
             - Returns: `Vector3D` - the scaled vector.

### 3. `NurbsCurve` Class
    - **Description**: Represents a NURBS curve in 3D space.
    - **Methods**:
        - `constructor(degree, controlPoints, knots, weights)`: Creates a new NURBS curve object.
            -   `degree` (Number): The degree of the curve.
            -   `controlPoints` (Array of Point3D): The control points of the curve.
            -   `knots` (Array of Number): The knot vector of the curve.
            -   `weights` (Array of Number): The weights of the control points.
        - `getDegree()`: Returns the degree of the curve.
        - `getControlPoints()`: Returns the control points of the curve.
        - `getKnots()`: Returns the knot vector of the curve.
        - `getWeights()`: Returns the weights of the control points.
        - `evaluate(u)`: Evaluates the curve at parameter `u`.
            -  `u` (Number): The parameter at which to evaluate the curve.
            -  Returns: `Point3D`: The point on the curve at parameter `u`.
        - `evaluateDerivative(u, order)`: Evaluates the derivatives of the curve at parameter `u` up to a given order
            -   `u` (Number): The parameter at which to evaluate the derivative.
            -   `order` (Number): The order of the derivative to evaluate.
             -  Returns: `Vector3D`: The derivative of the curve at parameter `u`.
        - `length()`: Computes the length of the NURBS curve.
             - Returns: `Number` - the arc length of the curve.
        - `tangeant(u)`: Returns the tangent vector at the parameter `u` of the curve.
              - `u` (Number): The parameter at which to evaluate the tangent.
              - Returns: `Vector3D` - the tangeant at point `u`.
        - `insertKnot(u, multiplicity)`: Inserts a knot at parameter u with a given multiplicity.
            -   `u` (Number): The parameter at which to insert the knot.
            -   `multiplicity` (Number): The number of times to insert the knot.
        - `removeKnot(u, multiplicity)`: Removes a knot at parameter u with a given multiplicity.
            -   `u` (Number): The parameter at which to remove the knot.
            -   `multiplicity` (Number): The number of times to remove the knot.
        -  `splitCurve(u)`: Splits the curve at parameter `u` into two curves
             -  `u` (Number): The parameter at which to split the curve.
             - Returns: `Array<NurbsCurve>` - an array of two nurbs curves.
         - `reverse()`: returns a new reversed NurbsCurve
             - Returns: `NurbsCurve` - the reversed curve.

### 4. `NurbsSurface` Class
    - **Description**: Represents a NURBS surface in 3D space.
    - **Methods**:
        - `constructor(degreeU, degreeV, controlPoints, knotsU, knotsV, weights)`: Creates a new NURBS surface object.
            -   `degreeU` (Number): The degree of the surface in the U direction.
            -   `degreeV` (Number): The degree of the surface in the V direction.
            -   `controlPoints` (Array of Array of Point3D): The control points of the surface.
            -   `knotsU` (Array of Number): The knot vector of the surface in the U direction.
            -   `knotsV` (Array of Number): The knot vector of the surface in the V direction.
            -   `weights` (Array of Array of Number): The weights of the control points.
        - `getDegreeU()`: Returns the degree of the surface in the U direction.
        - `getDegreeV()`: Returns the degree of the surface in the V direction.
        - `getControlPoints()`: Returns the control points of the surface.
        - `getKnotsU()`: Returns the knot vector of the surface in the U direction.
        - `getKnotsV()`: Returns the knot vector of the surface in the V direction.
        - `getWeights()`: Returns the weights of the control points.
        - `evaluate(u, v)`: Evaluates the surface at parameters `u` and `v`.
             -  `u` (Number): The parameter in the U direction.
             -  `v` (Number): The parameter in the V direction.
             -  Returns: `Point3D`: The point on the surface at parameters `u` and `v`.
         - `evaluatePartialDerivative(u, v, orderU, orderV)`: Evaluates the partial derivatives of the surface at parameters `u` and `v`
            - `u` (Number): The parameter in the U direction.
            - `v` (Number): The parameter in the V direction.
            - `orderU` (Number): The order of the derivative in the U direction.
            - `orderV` (Number): The order of the derivative in the V direction.
             -  Returns: `Vector3D`: The partial derivative of the surface at parameters `u` and `v`.
         -  `normal(u,v)`: Returns the normal vector at parameters `u` and `v`.
              - `u` (Number): The parameter in the U direction.
              - `v` (Number): The parameter in the V direction.
              - Returns: `Vector3D`: The normal vector at the given parameters.
         - `insertKnotU(u, multiplicity)`: Inserts a knot at parameter u in the U direction with a given multiplicity.
            -   `u` (Number): The parameter at which to insert the knot.
            -   `multiplicity` (Number): The number of times to insert the knot.
         - `insertKnotV(v, multiplicity)`: Inserts a knot at parameter v in the V direction with a given multiplicity.
            -   `v` (Number): The parameter at which to insert the knot.
            -   `multiplicity` (Number): The number of times to insert the knot.
         - `removeKnotU(u, multiplicity)`: Removes a knot at parameter u in the U direction with a given multiplicity.
            -   `u` (Number): The parameter at which to remove the knot.
            -   `multiplicity` (Number): The number of times to remove the knot.
         - `removeKnotV(v, multiplicity)`: Removes a knot at parameter v in the V direction with a given multiplicity.
            -   `v` (Number): The parameter at which to remove the knot.
            -   `multiplicity` (Number): The number of times to remove the knot.
         -  `splitSurfaceU(u)`: Splits the surface at parameter `u` in the u direction into two surfaces
             -  `u` (Number): The parameter at which to split the surface.
             - Returns: `Array<NurbsSurface>` - an array of two nurbs surfaces.
          - `splitSurfaceV(v)`: Splits the surface at parameter `v` in the v direction into two surfaces
             -  `v` (Number): The parameter at which to split the surface.
             - Returns: `Array<NurbsSurface>` - an array of two nurbs surfaces.

### 5. `Utils` Class
    - **Description**: Utility class providing functions for NURBS calculations.
    - **Methods**:
         - `basisFunction(i, degree, u, knots)`: Calculates the value of the i-th B-spline basis function at parameter `u`.
            -   `i` (Number): The index of the basis function.
            -   `degree` (Number): The degree of the B-spline.
            -   `u` (Number): The parameter value.
            -   `knots` (Array of Number): The knot vector.
            -   Returns: `Number` - the basis function value.
         - `basisFunctionDerivative(i, degree, u, knots, order)`: Calculates the derivative of the i-th B-spline basis function at parameter `u` of the given order.
            -   `i` (Number): The index of the basis function.
            -   `degree` (Number): The degree of the B-spline.
            -   `u` (Number): The parameter value.
            -    `knots` (Array of Number): The knot vector.
            -   `order` (Number): The order of the derivative.
            -    Returns: `Number` - the derivative of the basis function value.
         - `deBoor(degree, controlPoints, knots, weights, u)`: Evaluates the point on the curve at parameter `u` using the de Boor algorithm.
            -   `degree` (Number): The degree of the B-spline.
            -   `controlPoints` (Array of Point3D): The control points of the curve.
            -   `knots` (Array of Number): The knot vector of the curve.
            -   `weights` (Array of Number): The weights of the control points.
            -   `u` (Number): The parameter at which to evaluate the curve.
            -   Returns: `Point3D`: The point on the curve.
        - `deBoorDerivative(degree, controlPoints, knots, weights, u, order)`: Evaluates the derivatives of the curve using the de Boor algorithm.
            -   `degree` (Number): The degree of the B-spline.
            -   `controlPoints` (Array of Point3D): The control points of the curve.
            -   `knots` (Array of Number): The knot vector of the curve.
            -   `weights` (Array of Number): The weights of the control points.
            -   `u` (Number): The parameter at which to evaluate the curve.
            -   `order` (Number): The order of the derivative to evaluate.
            -   Returns: `Vector3D` - The derivative of the curve.
        - `findSpan(n, degree, u, knots)`: Finds the knot span index for parameter `u`.
             -   `n` (Number): The number of control points.
            -   `degree` (Number): The degree of the B-spline.
             -   `u` (Number): The parameter value.
            -   `knots` (Array of Number): The knot vector.
             -   Returns: `Number` - the knot span index.

## Additional Notes

-   **Error Handling:** Each method should include appropriate error handling and validation of inputs.
-   **Efficiency:** Consider performance implications when implementing algorithms, especially for high-degree curves/surfaces.
-   **Modularity:** Aim for a modular design, allowing for easy extension and modification of the library in the future.
-   **Documentation:** Clear and concise documentation is crucial for usability.
```