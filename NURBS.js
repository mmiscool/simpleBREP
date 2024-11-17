import * as THREE from 'three';
class NURBS_Curve {
    constructor(degree, controlPoints, knots) {
        this.degree = degree;
        this.controlPoints = controlPoints.map(cp => [
            ...cp,
            0
        ]);
        this.knots = knots;
        this.validateKnotVector();
    }
    derivative(t) {
        const degree = this.degree;
        const span = this.findSpan(degree, t, this.knots);
        const basisDeriv = this.basisFunctionDerivative(span, degree, t);
        const derivativePoint = new Array(this.controlPoints[0].length).fill(0);
        for (let i = 0; i < degree; i++) {
            const controlPoint = this.controlPoints[span - degree + i + 1];
            VectorUtils.addVectors(derivativePoint, VectorUtils.multiplyScalar(controlPoint, basisDeriv[i]));
        }
        return derivativePoint;
    }
    getPoints(stepSize) {
        if (stepSize <= 0 || stepSize > 1) {
            throw new Error('Step size must be a positive value between 0 and 1.');
        }
        const points = [];
        for (let t = 0; t <= 1; t += stepSize) {
            points.push(this.evaluate(t));
        }
        // To include the last point if step size does not land exactly on 1
        if (Math.abs(points[points.length - 1][0] - 1) > Number.EPSILON) {
            points.push(this.evaluate(1));
        }
        return points;
    }
    evaluate(t) {
        const point = new Array(this.controlPoints[0].length).fill(0);
        const degree = this.degree;
        const span = this.findSpan(degree, t, this.knots);
        for (let i = 0; i <= degree; i++) {
            const basis = this.basisFunction(span - degree + i, degree, t);
            const controlPoint = this.controlPoints[span - degree + i];
            VectorUtils.addVectors(point, VectorUtils.multiplyScalar(controlPoint, basis));
        }
        return point;
    }
    normal(t) {
        const derivative = this.derivative(t);
        const normalizedDerivative = VectorUtils.normalizeVector(derivative);
        const referenceVector = [
            1,
            0,
            0
        ];
        return this.crossProduct(normalizedDerivative, referenceVector);
    }
    insertKnot(knot) {
        const {knots, degree} = this;
        const newKnots = [];
        let insertIndex = knots.findIndex(k => k > knot);
        if (insertIndex === -1) {
            insertIndex = knots.length;
        }
        for (let i = 0; i < insertIndex; i++) {
            newKnots.push(knots[i]);
        }
        newKnots.push(knot);
        for (let i = insertIndex; i < knots.length; i++) {
            newKnots.push(knots[i]);
        }
        const newControlPoints = [...this.controlPoints];
        const multiplicity = this.knotMultiplicity(knot);
        const span = this.findSpan(degree, knot, newKnots);
        if (multiplicity < degree) {
            for (let j = 0; j < multiplicity; j++) {
                const basis = this.basisFunction(span, degree, knot);
                newControlPoints.splice(span - degree + 1, 0, this.buildControlPoint(span, knot, basis, newControlPoints));
            }
        }
        this.updateKnotVectorAndControlPoints(newKnots, newControlPoints);
    }
    removeKnot(knot) {
        const {knots} = this;
        const index = knots.indexOf(knot);
        if (index === -1) {
            throw new Error(`Knot ${ knot } not found in the knot vector.`);
        }
        const newKnots = [...knots];
        const multiplicity = this.knotMultiplicity(knot);
        if (multiplicity > 1) {
            for (let i = 0; i < multiplicity - 1; i++) {
                newKnots.splice(index, 1);
            }
        } else {
            newKnots.splice(index, 1);
            const span = this.findSpan(this.degree, knot, newKnots);
            for (let i = this.degree; i < this.controlPoints.length; i++) {
                const basis = this.basisFunction(span, this.degree, knot);
                this.controlPoints[i] = VectorUtils.subtractVectors(this.controlPoints[i], VectorUtils.multiplyScalar(this.controlPoints[i], basis));
            }
        }
        this.knots = newKnots;
    }
    knotMultiplicity(knot) {
        return this.knots.filter(currentKnot => currentKnot === knot).length;
    }
    degreeElevate() {
        const newDegree = this.degree + 1;
        const newControlPoints = [];
        const knotVector = this.knots.slice(this.degree, this.knots.length - this.degree);
        const newKnots = [
            ...knotVector,
            ...Array(newDegree + 1).fill(1)
        ];
        const n = this.controlPoints.length - 1;
        for (let i = 0; i < n; i++) {
            const newPoint = new Array(this.controlPoints[0].length).fill(0);
            for (let j = 0; j <= this.degree; j++) {
                const span = this.findSpan(this.degree, knotVector[i], this.knots);
                const basis = this.basisFunction(span, this.degree, knotVector[i]);
                newPoint = VectorUtils.addVectors(newPoint, VectorUtils.multiplyScalar(this.controlPoints[span - this.degree + j], basis));
            }
            newControlPoints.push(newPoint);
        }
        this.degree = newDegree;
        this.controlPoints = newControlPoints;
        this.knots = newKnots;
    }
    basisFunction(i, k, t) {
        if (k === 0) {
            return t >= this.knots[i] && t < this.knots[i + 1] ? 1 : 0;
        }
        const denominator1 = this.knots[i + k] - this.knots[i];
        const denominator2 = this.knots[i + k + 1] - this.knots[i + 1];
        const left = this.basisFunction(i, k - 1, t);
        const right = this.basisFunction(i + 1, k - 1, t);
        const term1 = denominator1 !== 0 ? (t - this.knots[i]) / denominator1 * left : 0;
        const term2 = denominator2 !== 0 ? (this.knots[i + k + 1] - t) / denominator2 * right : 0;
        return term1 + term2;
    }
    length(stepSize = 0.01) {
        if (stepSize <= 0) {
            throw new Error('Step size must be a positive value.');
        }
        let totalLength = 0;
        let previousPoint = this.evaluate(0);
        for (let t = stepSize; t <= 1; t += stepSize) {
            const currentPoint = this.evaluate(t);
            totalLength += VectorUtils.distanceBetweenPoints(currentPoint, previousPoint);
            previousPoint = currentPoint;
        }
        return totalLength;
    }
    approximate(points) {
        const degree = this.degree;
        const n = points.length;
        const newControlPoints = Array.from({ length: degree + 1 }, () => new Array(points[0].length).fill(0));
        const stepSize = 1 / (n - 1);
        for (let i = 0; i <= degree; i++) {
            for (let j = 0; j < n; j++) {
                const t = j * stepSize;
                const basis = this.basisFunction(i, degree, t);
                newControlPoints[i] = VectorUtils.addVectors(newControlPoints[i], VectorUtils.multiplyScalar(points[j], basis));
            }
        }
        this.controlPoints = newControlPoints;
    }
    normalizeParameterization() {
        const n = this.controlPoints.length - 1;
        const newKnots = [];
        const newControlPoints = [...this.controlPoints];
        const totalLength = this.length();
        if (totalLength === 0) {
            throw new Error('Total length of the curve is zero. Cannot normalize.');
        }
        const accumulatedLength = new Array(n + 1).fill(0);
        for (let i = 0; i < n; i++) {
            const p0 = this.evaluate(i / n);
            const p1 = this.evaluate((i + 1) / n);
            accumulatedLength[i + 1] = accumulatedLength[i] + VectorUtils.distanceBetweenPoints(p0, p1);
        }
        for (let i = 0; i <= n; i++) {
            const normalizedParameter = accumulatedLength[i] / totalLength;
            newKnots.push(i === 0 || i === n ? normalizedParameter : normalizedParameter);
        }
        this.knots = newKnots;
        this.controlPoints = newControlPoints;
    }
    findSpan(degree, t, knots) {
        const n = knots.length - degree - 1;
        if (t >= knots[n + 1])
            return n;
        if (t <= knots[degree])
            return degree;
        let low = degree;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (t < knots[mid] || t >= knots[mid + 1]) {
            if (t < knots[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    getKnotVector() {
        return this.knots.slice();
    }
    static fromControlPoints(degree, controlPoints) {
        const n = controlPoints.length;
        const knots = [];
        for (let i = 0; i <= degree; i++) {
            knots.push(0);
        }
        for (let i = 0; i < n - degree; i++) {
            knots.push(i / (n - degree));
        }
        for (let i = 0; i <= degree; i++) {
            knots.push(1);
        }
        return new NURBS_Curve(degree, controlPoints, knots);
    }
    selfIntersect() {
        const intersections = [];
        const stepSize = 0.01;
        const evaluatedPoints = [];
        const evaluatedTValues = [];
        for (let t = 0; t <= 1; t += stepSize) {
            evaluatedPoints.push(this.evaluate(t));
            evaluatedTValues.push(t);
        }
        const n = evaluatedPoints.length;
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                if (this.isClose(evaluatedPoints[i], evaluatedPoints[j])) {
                    intersections.push({
                        point: evaluatedPoints[i],
                        t1: evaluatedTValues[i],
                        t2: evaluatedTValues[j]
                    });
                }
            }
        }
        return intersections;
    }
    split(t) {
        const span = this.findSpan(this.degree, t, this.knots);
        const leftNewControlPoints = [];
        const rightNewControlPoints = [];
        for (let i = 0; i <= this.degree; i++) {
            leftNewControlPoints.push(this.controlPoints[span - this.degree + i]);
            rightNewControlPoints.push(this.controlPoints[span - this.degree + i]);
        }
        const leftKnots = [
            ...this.knots.slice(0, span - this.degree + 1),
            t
        ];
        const rightKnots = [
            t,
            ...this.knots.slice(span + 1)
        ];
        for (let i = 0; i < this.degree; i++) {
            const newPointLeft = new Array(this.controlPoints[0].length).fill(0);
            const newPointRight = new Array(this.controlPoints[0].length).fill(0);
            for (let j = 0; j <= this.degree; j++) {
                const leftBasis = this.basisFunction(span - this.degree + j, this.degree, t);
                const rightBasis = this.basisFunction(span - this.degree + j + 1, this.degree, t);
                VectorUtils.addVectors(newPointLeft, VectorUtils.multiplyScalar(this.controlPoints[span - this.degree + j], leftBasis));
                VectorUtils.addVectors(newPointRight, VectorUtils.multiplyScalar(this.controlPoints[span - this.degree + j + 1], rightBasis));
            }
            leftNewControlPoints.push(newPointLeft);
            rightNewControlPoints.push(newPointRight);
        }
        const leftCurve = new NURBS_Curve(this.degree, leftNewControlPoints, leftKnots);
        const rightCurve = new NURBS_Curve(this.degree, rightNewControlPoints, rightKnots);
        return {
            leftCurve,
            rightCurve
        };
    }
    isPointOnCurve(point, tolerance = 0.000001) {
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const curvePoint = this.evaluate(t);
            if (this.isClose(curvePoint, point, tolerance)) {
                return true;
            }
        }
        return false;
    }
    closestPoint(point) {
        let minDistance = Infinity;
        let closestPoint = null;
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const curvePoint = this.evaluate(t);
            const distance = VectorUtils.distanceBetweenPoints(curvePoint, point);
            if (distance < minDistance) {
                minDistance = distance;
                closestPoint = curvePoint;
            }
        }
        return closestPoint;
    }
    getCurveInfo() {
        return {
            degree: this.degree,
            controlPointsCount: this.controlPoints.length,
            knotsCount: this.knots.length
        };
    }
    intersectsWith(otherCurve) {
        const intersections = [];
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const pointA = this.evaluate(t);
            for (let s = 0; s <= 1; s += stepSize) {
                const pointB = otherCurve.evaluate(s);
                if (this.isClose(pointA, pointB)) {
                    intersections.push({
                        point: pointA,
                        paramA: t,
                        paramB: s
                    });
                }
            }
        }
        return intersections;
    }
    uniformSubdivision(n) {
        const points = [];
        for (let i = 0; i <= n; i++) {
            const t = i / n;
            points.push(this.evaluate(t));
        }
        return points;
    }
    convexHull() {
        const hull = [];
        const len = this.controlPoints.length;
        for (let i = 0; i < len; i++) {
            hull.push(this.controlPoints[i]);
        }
        return hull;
    }
    centroid() {
        const n = this.controlPoints.length;
        const sum = this.controlPoints.reduce((acc, point) => VectorUtils.addVectors(acc, point), new Array(this.controlPoints[0].length).fill(0));
        return sum.map(coord => coord / n);
    }
    isClosed() {
        return this.controlPoints[0] === this.controlPoints[this.controlPoints.length - 1];
    }
    arcLength(stepSize = 0.01) {
        if (stepSize <= 0) {
            throw new Error('Step size must be a positive value.');
        }
        let totalLength = 0;
        let previousPoint = this.evaluate(0);
        for (let t = stepSize; t <= 1; t += stepSize) {
            const currentPoint = this.evaluate(t);
            totalLength += VectorUtils.distanceBetweenPoints(currentPoint, previousPoint);
            previousPoint = currentPoint;
        }
        return totalLength;
    }
    basisFunctionDerivative(i, k, t) {
        const left = new Array(k).fill(0);
        const right = new Array(k).fill(0);
        const result = new Array(k - 1).fill(0);
        for (let j = 0; j < k; j++) {
            left[j] = t - this.knots[i - j];
            right[j] = this.knots[i + 1 + j] - t;
        }
        for (let j = 0; j < k - 1; j++) {
            const denominator = right[j + 1] + left[j];
            if (denominator !== 0) {
                const temp = (k - 1) / denominator;
                result[j] = temp * (right[j + 1] * this.basisFunctionDerivative(i, k - 1, t) + left[j] * (j > 0 ? result[j - 1] : 0));
            }
        }
        return result;
    }
    crossProduct(v1, v2) {
        if (v1.length !== 3 || v2.length !== 3) {
            throw new Error('Both input vectors must be 3D.');
        }
        return [
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        ];
    }
    isClose(point1, point2, tolerance = 0.000001) {
        return point1.every((val, index) => Math.abs(val - point2[index]) < tolerance);
    }
    area(stepSize = 0.01) {
        if (stepSize <= 0) {
            throw new Error('Step size must be a positive value.');
        }
        let totalArea = 0;
        const evaluatedPoints = this.getPoints(stepSize);
        for (let i = 0; i < evaluatedPoints.length - 1; i++) {
            const p0 = evaluatedPoints[i];
            const p1 = evaluatedPoints[i + 1];
            const p2 = [
                0,
                0,
                0
            ];
            totalArea += this.triangleArea(p0, p1, p2);
        }
        return totalArea;
    }
    trimCurve(curve) {
        const trimmedControlPoints = [];
        const trimmedKnots = [...this.knots];
        const controlPointSet = new Set();
        const trimPoints = [];
        const stepSize = (curve.knots.length - 1) / 100;
        for (let t = curve.knots[curve.degree]; t <= curve.knots[curve.knots.length - curve.degree - 1]; t += stepSize) {
            const trimPoint = curve.evaluate(t);
            trimPoints.push(trimPoint);
        }
        for (const controlPoint of this.controlPoints) {
            if (this.isControlPointInsideTrimCurve(controlPoint, trimPoints)) {
                controlPointSet.add(controlPoint);
            }
        }
        this.controlPoints = this.controlPoints.filter(controlPoint => controlPointSet.has(controlPoint));
        this.controlPoints.forEach((cp, index) => {
            this.controlPoints[index] = VectorUtils.addVectors(cp, [
                0,
                0,
                0
            ]);
        });
        this.knots = trimmedKnots;
    }
    fitToPoints(points) {
        const degree = this.degree;
        const n = points.length;
        const newControlPoints = Array.from({ length: n }, () => new Array(points[0].length).fill(0));
        const stepSize = 1 / (n - 1);
        const knotVector = this.parameterizeKnots(n, degree);
        for (let i = 0; i < n; i++) {
            const t = i * stepSize;
            for (let k = 0; k <= degree; k++) {
                const basisValue = this.basisFunction(k, degree, t);
                for (let idx = 0; idx < newControlPoints[k].length; idx++) {
                    newControlPoints[k][idx] += basisValue * points[i][idx];
                }
            }
        }
        this.controlPoints = newControlPoints;
        this.knots = knotVector;
    }
    degreeElevateU() {
        const newDegreeU = this.degreeU + 1;
        const newControlPoints = [];
        const knotVectorU = this.knotsU.slice(this.degreeU, this.knotsU.length - this.degreeU);
        for (let i = 0; i <= this.controlPoints.length - 1; i++) {
            const newPoint = new Array(this.controlPoints[0].length).fill(0);
            for (let j = 0; j <= this.degreeU; j++) {
                const span = this.findSpanU(knotVectorU[i], this.degreeU);
                const basis = this.basisFunctionU(span, this.degreeU, knotVectorU[i]);
                for (let k = 0; k < newPoint.length; k++) {
                    newPoint[k] += basis * this.controlPoints[span - this.degreeU + j][k];
                }
            }
            newControlPoints.push(newPoint);
        }
        this.degreeU = newDegreeU;
        this.controlPoints = newControlPoints;
        this.parameterize();
    }
    degreeElevateV() {
        const newDegreeV = this.degreeV + 1;
        const newControlPoints = [...this.controlPoints];
        // Update knot vector and control points logic
        this.degreeV = newDegreeV;
        this.controlPoints = newControlPoints;
    }
    parameterize() {
        const n = this.controlPoints.length;
        const newKnots = [];
        for (let i = 0; i <= this.degree; i++) {
            newKnots.push(0);
        }
        for (let i = 0; i < n - this.degree; i++) {
            newKnots.push(i / (n - this.degree));
        }
        for (let i = 0; i <= this.degree; i++) {
            newKnots.push(1);
        }
        this.knots = newKnots;
    }
    surfaceArea(stepSizeU, stepSizeV) {
        let totalArea = 0;
        const points = this.getPoints(stepSizeU, stepSizeV);
        for (let i = 0; i < points.length - 1; i++) {
            const p0 = points[i];
            const p1 = points[i + 1];
            const p2 = points[(i + 1) % points.length];
            totalArea += this.triangleArea(p0, p1, p2);
        }
        return totalArea;
    }
    parameterizeKnots(n, degree) {
        const knots = [];
        for (let i = 0; i <= degree; i++) {
            knots.push(0);
        }
        for (let i = 0; i < n - degree; i++) {
            knots.push(i / (n - degree));
        }
        for (let i = 0; i <= degree; i++) {
            knots.push(1);
        }
        return knots;
    }
    computeFirstDerivatives(stepSize = 0.01) {
        const derivatives = [];
        for (let t = 0; t <= 1; t += stepSize) {
            const derivative = this.derivative(t);
            derivatives.push(derivative);
        }
        return derivatives;
    }
    isControlPointInsideTrimCurve(controlPoint, trimPoints) {
        const [x, y] = controlPoint;
        const minX = Math.min(...trimPoints.map(p => p[0]));
        const maxX = Math.max(...trimPoints.map(p => p[0]));
        const minY = Math.min(...trimPoints.map(p => p[1]));
        const maxY = Math.max(...trimPoints.map(p => p[1]));
        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }
    triangleArea(p0, p1, p2) {
        const a = VectorUtils.distanceBetweenPoints(p0, p1);
        const b = VectorUtils.distanceBetweenPoints(p1, p2);
        const c = VectorUtils.distanceBetweenPoints(p2, p0);
        const s = (a + b + c) / 2;
        return Math.sqrt(s * (s - a) * (s - b) * (s - c));
    }
    distanceBetweenPoints(p1, p2) {
        return Math.sqrt(VectorUtils.subtractVectors(p1, p2).reduce((sum, val) => sum + val ** 2, 0));
    }
    reverse() {
        this.controlPoints.reverse();
        this.knots.reverse();
    }
    getControlPoint(index) {
        if (index < 0 || index >= this.controlPoints.length) {
            throw new Error(`Control point index ${ index } is out of bounds.`);
        }
        return this.controlPoints[index];
    }
    updateKnotVectorAndControlPoints(newKnots, newControlPoints) {
        if (!Array.isArray(newKnots) || !Array.isArray(newControlPoints)) {
            throw new Error('New knots and control points must be arrays.');
        }
        if (newKnots.length < this.degree + 2) {
            throw new Error(`New knot vector must contain at least ${ this.degree + 2 } elements.`);
        }
        if (newControlPoints.length !== this.controlPoints.length) {
            throw new Error('New control points array must match the length of the current control points.');
        }
        this.knots = newKnots;
        this.controlPoints = newControlPoints;
    }
    calculateCurvature(pointA, pointB, pointC) {
        const tangentAB = this.calculateTangent(pointA, pointB);
        const tangentAC = this.calculateTangent(pointA, pointC);
        const curvature = this.crossProduct(tangentAB, tangentAC);
        const magnitude = Math.sqrt(curvature.reduce((sum, val) => sum + val * val, 0));
        if (magnitude === 0) {
            throw new Error('Curvature at the specified points results in a zero vector.');
        }
        return curvature.map(val => val / magnitude);
    }
    blendControlPoints(pointA, pointB, knot, span, degree) {
        const newPoint = new Array(pointA.length).fill(0);
        const basisA = this.basisFunction(span, degree, knot);
        const basisB = this.basisFunction(span - 1, degree, knot);
        for (let i = 0; i < pointA.length; i++) {
            newPoint[i] = basisA * pointA[i] + basisB * pointB[i];
        }
        return newPoint;
    }
    validateKnotVector() {
        const n = this.controlPoints.length;
        const knotMultiplicity = this.knots.reduce((acc, knot) => {
            acc[knot] = (acc[knot] || 0) + 1;
            return acc;
        }, {});
        const minKnot = this.knots[0];
        const maxKnot = this.knots[this.knots.length - 1];
        if (this.knots.length < this.degree + n + 1) {
            throw new Error('Knot vector must have at least degree + number of control points + 1 elements.');
        }
        if (minKnot !== 0 || maxKnot !== 1) {
            throw new Error('Knot vector must start with 0 and end with 1.');
        }
        for (const knot in knotMultiplicity) {
            if (knotMultiplicity[knot] > this.degree + 1) {
                throw new Error(`Knot ${ knot } has a multiplicity greater than degree + 1.`);
            }
        }
    }
    reorderKnotVector(newKnots) {
        this.knots = newKnots;
        this.validateKnotVector();
    }
    evaluateSurfaceDerivative(u, v, direction) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const derivativePoint = new Array(this.controlPoints[0][0].length).fill(0);
        if (direction === 'u') {
            for (let i = 0; i <= this.degreeU; i++) {
                for (let j = 0; j <= this.degreeV; j++) {
                    const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                    const basisDerivativeU = this.basisFunctionU(spanU, this.degreeU - 1, u);
                    for (let k = 0; k < derivativePoint.length; k++) {
                        derivativePoint[k] += basisDerivativeU * basisV[j] * controlPoint[k];
                    }
                }
            }
        } else if (direction === 'v') {
            for (let i = 0; i <= this.degreeU; i++) {
                for (let j = 0; j <= this.degreeV; j++) {
                    const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                    const basisDerivativeV = this.basisFunctionV(spanV, this.degreeV - 1, v);
                    for (let k = 0; k < derivativePoint.length; k++) {
                        derivativePoint[k] += basisU[i] * basisDerivativeV * controlPoint[k];
                    }
                }
            }
        } else {
            throw new Error('Invalid direction specified. Use \'u\' or \'v\'.');
        }
        return derivativePoint;
    }
    basisFunctionU(i, k, u) {
        if (k === 0) {
            return u >= this.knotsU[i] && u < this.knotsU[i + 1] ? 1 : 0;
        }
        const left = this.basisFunctionU(i, k - 1, u);
        const right = this.basisFunctionU(i + 1, k - 1, u);
        const denominator1 = this.knotsU[i + k] - this.knotsU[i];
        const denominator2 = this.knotsU[i + k + 1] - this.knotsU[i + 1];
        const term1 = denominator1 !== 0 ? (u - this.knotsU[i]) / denominator1 * left : 0;
        const term2 = denominator2 !== 0 ? (this.knotsU[i + k + 1] - u) / denominator2 * right : 0;
        return term1 + term2;
    }
    isPointInsideTrimCurve(point, trimPoints) {
        const [x, y] = point;
        const minX = Math.min(...trimPoints.map(p => p[0]));
        const maxX = Math.max(...trimPoints.map(p => p[0]));
        const minY = Math.min(...trimPoints.map(p => p[1]));
        const maxY = Math.max(...trimPoints.map(p => p[1]));
        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }
    generateControlPoints() {
        const numPoints = 100;
        for (let i = 0; i < numPoints; i++) {
            const angle = i / numPoints * 2 * Math.PI;
            const x = this.center[0] + this.radius * Math.cos(angle);
            const y = this.center[1] + this.radius * Math.sin(angle);
            this.controlPoints.push(VectorUtils.addVectors([
                x,
                y,
                this.center[2]
            ], [
                0,
                0,
                0
            ]));
        }
    }
    secondDerivative(t) {
        const degree = this.degree;
        const span = this.findSpan(degree, t, this.knots);
        const basisDeriv = this.basisFunctionDerivative(span, degree, t);
        let secondDerivPoint = new Array(this.controlPoints[0].length).fill(0);
        for (let i = 0; i < degree - 1; i++) {
            const controlPoint = this.controlPoints[span - degree + i + 1];
            secondDerivPoint = VectorUtils.addVectors(secondDerivPoint, VectorUtils.multiplyScalar(controlPoint, basisDeriv[i]));
        }
        return secondDerivPoint;
    }
    tangent(t) {
        const derivativePoint = this.derivative(t);
        return VectorUtils.normalizeVector(derivativePoint);
    }
    /* *
     * Calculates the curvature of the NURBS curve at a specific parameter t.
     * @param {number} t - The parameter value along the curve.
     * @returns {Array<number>} - A vector representing the curvature.*/
    curvature(t) {
        const firstDerivative = this.derivative(t);
        const secondDerivative = this.secondDerivative(t);
        const curvatureVector = this.crossProduct(firstDerivative, secondDerivative);
        const magnitude = Math.sqrt(curvatureVector.reduce((sum, val) => sum + val * val, 0));
        if (magnitude === 0) {
            throw new Error('Curvature results in a zero vector. Cannot compute curvature.');
        }
        return curvatureVector.map(val => val / magnitude);
    }
    getTangent(t) {
        const derivative = this.derivative(t);
        const length = Math.sqrt(derivative.reduce((acc, val) => acc + val * val, 0));
        if (length === 0) {
            throw new Error('The derivative at the parameter t results in a zero vector. Cannot compute tangent.');
        }
        return derivative.map(val => val / length);
    }
    addVectors(v1, v2) {
        return VectorUtils.addVectors(v1, v2);
    }
    // Utility method for vector subtraction
    static subtractVectors(v1, v2) {
        return v1.map((val, index) => val - v2[index]);
    }
    // Utility method for vector normalization
    static normalizeVector(v) {
        const length = Math.sqrt(v.reduce((sum, val) => sum + val ** 2, 0));
        if (length === 0) {
            throw new Error('Cannot normalize a zero-length vector.');
        }
        return v.map(val => val / length);
    }
    nearestPoint(point) {
        let minDistance = Infinity;
        let closestPoint = null;
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const curvePoint = this.evaluate(t);
            const distance = VectorUtils.distanceBetweenPoints(curvePoint, point);
            if (distance < minDistance) {
                minDistance = distance;
                closestPoint = curvePoint;
            }
        }
        return closestPoint;
    }
    intersectionsWith(otherCurve) {
        const intersections = [];
        const stepSize = 0.01;
        const evaluatedPointsA = [];
        const evaluatedPointsB = [];
        for (let t = 0; t <= 1; t += stepSize) {
            evaluatedPointsA.push(this.evaluate(t));
        }
        for (let t = 0; t <= 1; t += stepSize) {
            evaluatedPointsB.push(otherCurve.evaluate(t));
        }
        for (let i = 0; i < evaluatedPointsA.length; i++) {
            for (let j = 0; j < evaluatedPointsB.length; j++) {
                if (this.isClose(evaluatedPointsA[i], evaluatedPointsB[j])) {
                    intersections.push({
                        point: evaluatedPointsA[i],
                        t1: i * stepSize,
                        t2: j * stepSize
                    });
                }
            }
        }
        return intersections;
    }
    normalize(t) {
        const derivativePoint = this.derivative(t);
        return VectorUtils.normalizeVector(derivativePoint);
    }
    buildControlPoint(span, t, basis, newControlPoints) {
        const controlPoint = new Array(this.controlPoints[0].length).fill(0);
        for (let i = 0; i < this.degree + 1; i++) {
            const cp = this.controlPoints[span - this.degree + i];
            for (let j = 0; j < controlPoint.length; j++) {
                controlPoint[j] += cp[j] * basis[i];
            }
        }
        return controlPoint;
    }
    generateControlPoint(span, t, basis, newControlPoints) {
        const controlPoint = this.buildControlPoint(span, t, basis, newControlPoints);
        newControlPoints.push(controlPoint);
    }
    toBufferGeometry(stepSize = 0.01) {
        const points = this.getPoints(stepSize);
        const geometry = new THREE.BufferGeometry();
        const vertices = new Float32Array(points.flat());
        geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
        return geometry;
    }
    closestIntersection(otherCurve) {
        const intersections = this.intersectsWith(otherCurve);
        if (intersections.length === 0)
            return null;
        let minDistance = Infinity;
        let closestIntersection = null;
        intersections.forEach(({point, paramA, paramB}) => {
            const distance = VectorUtils.distanceBetweenPoints(point, this.evaluate(paramA));
            if (distance < minDistance) {
                minDistance = distance;
                closestIntersection = {
                    point,
                    paramA,
                    paramB
                };
            }
        });
        return closestIntersection;
    }
    basisFunctionV(j, m, v) {
        if (m === 0) {
            return v >= this.knotsV[j] && v < this.knotsV[j + 1] ? 1 : 0;
        }
        const left = this.basisFunctionV(j, m - 1, v);
        const right = this.basisFunctionV(j + 1, m - 1, v);
        const denominator1 = this.knotsV[j + m] - this.knotsV[j];
        const denominator2 = this.knotsV[j + m + 1] - this.knotsV[j + 1];
        const term1 = denominator1 !== 0 ? (v - this.knotsV[j]) / denominator1 * left : 0;
        const term2 = denominator2 !== 0 ? (this.knotsV[j + m + 1] - v) / denominator2 * right : 0;
        return term1 + term2;
    }
    findSpanU(u) {
        const n = this.knotsU.length - this.degreeU - 1;
        if (u >= this.knotsU[n + 1])
            return n;
        if (u <= this.knotsU[this.degreeU])
            return this.degreeU;
        let low = this.degreeU;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (u < this.knotsU[mid] || u >= this.knotsU[mid + 1]) {
            if (u < this.knotsU[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    findSpanV(v) {
        const n = this.knotsV.length - this.degreeV - 1;
        if (v >= this.knotsV[n + 1])
            return n;
        if (v <= this.knotsV[this.degreeV])
            return this.degreeV;
        let low = this.degreeV;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (v < this.knotsV[mid] || v >= this.knotsV[mid + 1]) {
            if (v < this.knotsV[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    evaluateSurface(u, v) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const point = new Array(this.controlPoints[0].length).fill(0);
        for (let i = 0; i <= this.degreeU; i++) {
            for (let j = 0; j <= this.degreeV; j++) {
                const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                const weightedPoint = VectorUtils.multiplyScalar(controlPoint, basisU[i] * basisV[j]);
                point.forEach((val, k) => point[k] += weightedPoint[k]);
            }
        }
        return point;
    }
    calculateBoundingBox() {
        const points = this.getPoints(0.01);
        const xValues = points.map(point => point[0]);
        const yValues = points.map(point => point[1]);
        const zValues = points.map(point => point[2]);
        const boundingBox = {
            min: [
                Math.min(...xValues),
                Math.min(...yValues),
                Math.min(...zValues)
            ],
            max: [
                Math.max(...xValues),
                Math.max(...yValues),
                Math.max(...zValues)
            ]
        };
        return boundingBox;
    }
    evaluateAtMultipleParameters(tValues) {
        return tValues.map(t => this.evaluate(t));
    }
    getControlPoints() {
        return this.controlPoints.slice();
    }
    getDegree() {
        return this.degree;
    }
    calculateCentroid() {
        const numControlPoints = this.controlPoints.length;
        const sum = this.controlPoints.reduce((acc, cp) => {
            return VectorUtils.addVectors(acc, cp);
        }, [
            0,
            0,
            0
        ]);
        return sum.map(coord => coord / numControlPoints);
    }
}
class NURBS_Surface {
    constructor(degreeU, degreeV, controlPoints, knotsU, knotsV) {
        this.degreeU = degreeU;
        this.degreeV = degreeV;
        this.controlPoints = controlPoints;
        this.knotsU = knotsU;
        this.knotsV = knotsV;
        this.validateKnotVectors();
    }
    evaluate(u, v) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const point = new Array(this.controlPoints[0][0].length).fill(0);
        for (let i = 0; i <= this.degreeU; i++) {
            for (let j = 0; j <= this.degreeV; j++) {
                const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                const weightedPoint = VectorUtils.multiplyScalar(controlPoint, basisU[i] * basisV[j]);
                point.forEach((val, k) => point[k] += weightedPoint[k]);
            }
        }
        return point;
    }
    derivative(u, v) {
        const derivativeU = this.evaluateSurfaceDerivative(u, v, 'u');
        const derivativeV = this.evaluateSurfaceDerivative(u, v, 'v');
        return {
            derivativeU,
            derivativeV
        };
    }
    getPoints(stepSizeU, stepSizeV) {
        const points = [];
        for (let u = 0; u <= 1; u += stepSizeU) {
            for (let v = 0; v <= 1; v += stepSizeV) {
                points.push(this.evaluate(u, v));
            }
        }
        return points;
    }
    normal(u, v) {
        const {derivativeU, derivativeV} = this.derivative(u, v);
        const normal = this.crossProduct(derivativeU, derivativeV);
        const length = Math.sqrt(normal.reduce((acc, val) => acc + val * val, 0));
        if (length === 0) {
            throw new Error('The computed normal is a zero vector at parameterization (u, v).');
        }
        return normal.map(val => val / length);
    }
    insertKnotU(knot) {
        const {knotsU} = this;
        const newKnotsU = [];
        let insertIndex = knotsU.findIndex(k => k > knot);
        if (insertIndex === -1) {
            insertIndex = knotsU.length;
        }
        for (let i = 0; i < insertIndex; i++) {
            newKnotsU.push(knotsU[i]);
        }
        newKnotsU.push(knot);
        for (let i = insertIndex; i < knotsU.length; i++) {
            newKnotsU.push(knotsU[i]);
        }
        this.knotsU = newKnotsU;
    }
    insertKnotV(knot) {
        const {knotsV} = this;
        const newKnotsV = [];
        let insertIndex = knotsV.findIndex(k => k > knot);
        if (insertIndex === -1) {
            insertIndex = knotsV.length;
        }
        for (let i = 0; i < insertIndex; i++) {
            newKnotsV.push(knotsV[i]);
        }
        newKnotsV.push(knot);
        for (let i = insertIndex; i < knotsV.length; i++) {
            newKnotsV.push(knotsV[i]);
        }
        this.knotsV = newKnotsV;
    }
    removeKnotU(knot) {
        const index = this.knotsU.indexOf(knot);
        if (index === -1) {
            throw new Error(`Knot ${ knot } not found in the knot vector U.`);
        }
        const newKnotsU = [...this.knotsU];
        const multiplicity = this.knotMultiplicity(knot, this.knotsU);
        if (multiplicity > 1) {
            for (let i = 0; i < multiplicity - 1; i++) {
                newKnotsU.splice(index, 1);
            }
        } else {
            newKnotsU.splice(index, 1);
            const span = this.findSpanU(knot);
            // Update control points if necessary
            for (let i = this.degreeU; i < this.controlPoints.length; i++) {
                const basis = this.basisFunctionU(span, this.degreeU, knot);
                for (let j = 0; j < this.controlPoints[i].length; j++) {
                    this.controlPoints[i][j] -= basis * this.controlPoints[i][j];
                }
            }
        }
        this.knotsU = newKnotsU;
    }
    removeKnotV(knot) {
        const index = this.knotsV.indexOf(knot);
        if (index === -1) {
            throw new Error(`Knot ${ knot } not found in the knot vector V.`);
        }
        const newKnotsV = [...this.knotsV];
        const multiplicity = this.knotMultiplicity(knot, this.knotsV);
        if (multiplicity > 1) {
            for (let i = 0; i < multiplicity - 1; i++) {
                newKnotsV.splice(index, 1);
            }
        } else {
            newKnotsV.splice(index, 1);
            const span = this.findSpanV(knot);
            // Update control points if necessary
            for (let i = this.degreeV; i < this.controlPoints[0].length; i++) {
                const basis = this.basisFunctionV(span, this.degreeV, knot);
                for (let j = 0; j < this.controlPoints.length; j++) {
                    this.controlPoints[j][i] -= basis * this.controlPoints[j][i];
                }
            }
        }
        this.knotsV = newKnotsV;
    }
    setControlPoint(i, j, point) {
        if (i < 0 || i >= this.controlPoints.length) {
            throw new Error(`Control point index i (${ i }) is out of bounds.`);
        }
        if (j < 0 || j >= this.controlPoints[i].length) {
            throw new Error(`Control point index j (${ j }) is out of bounds.`);
        }
        if (point.length !== this.controlPoints[0][0].length) {
            throw new Error(`Control point must have ${ this.controlPoints[0][0].length } dimensions.`);
        }
        this.controlPoints[i][j] = point;
    }
    getControlPoint(i, j) {
        if (i < 0 || i >= this.controlPoints.length) {
            throw new Error(`Control point index i (${ i }) is out of bounds.`);
        }
        if (j < 0 || j >= this.controlPoints[i].length) {
            throw new Error(`Control point index j (${ j }) is out of bounds.`);
        }
        return this.controlPoints[i][j];
    }
    degreeElevateU() {
        const newDegreeU = this.degreeU + 1;
        const newControlPoints = [];
        const knotVectorU = this.knotsU.slice(this.degreeU, this.knotsU.length - this.degreeU);
        for (let i = 0; i <= this.controlPoints.length - 1; i++) {
            const newPoint = new Array(this.controlPoints[0][0].length).fill(0);
            for (let j = 0; j <= this.degreeU; j++) {
                const span = this.findSpanU(knotVectorU[i], this.degreeU);
                const basis = this.basisFunctionU(span, this.degreeU, knotVectorU[i]);
                for (let k = 0; k < newPoint.length; k++) {
                    newPoint[k] += basis * this.controlPoints[span - this.degreeU + j][k];
                }
            }
            newControlPoints.push(newPoint);
        }
        this.degreeU = newDegreeU;
        this.controlPoints = newControlPoints;
        this.parameterizeSurface();
    }
    degreeElevateV() {
        const newDegreeV = this.degreeV + 1;
        const newControlPoints = [];
        const knotVectorV = this.knotsV.slice(this.degreeV, this.knotsV.length - this.degreeV);
        for (let i = 0; i < this.controlPoints.length; i++) {
            const newPoint = new Array(this.controlPoints[0][0].length).fill(0);
            for (let j = 0; j <= this.degreeV; j++) {
                const span = this.findSpanV(knotVectorV[i], this.degreeV);
                const basis = this.basisFunctionV(span, this.degreeV, knotVectorV[i]);
                for (let k = 0; k < newPoint.length; k++) {
                    newPoint[k] += basis * this.controlPoints[span - this.degreeV + j][k];
                }
            }
            newControlPoints.push(newPoint);
        }
        this.degreeV = newDegreeV;
        this.controlPoints = newControlPoints;
        this.parameterizeSurface();
    }
    basisFunctionU(i, k, u) {
        if (k === 0) {
            return u >= this.knotsU[i] && u < this.knotsU[i + 1] ? 1 : 0;
        }
        const left = this.basisFunctionU(i, k - 1, u);
        const right = this.basisFunctionU(i + 1, k - 1, u);
        const denominator1 = this.knotsU[i + k] - this.knotsU[i];
        const denominator2 = this.knotsU[i + k + 1] - this.knotsU[i + 1];
        const term1 = denominator1 !== 0 ? (u - this.knotsU[i]) / denominator1 * left : 0;
        const term2 = denominator2 !== 0 ? (this.knotsU[i + k + 1] - u) / denominator2 * right : 0;
        return term1 + term2;
    }
    basisFunctionV(j, m, v) {
        if (m === 0) {
            return v >= this.knotsV[j] && v < this.knotsV[j + 1] ? 1 : 0;
        }
        const left = this.basisFunctionV(j, m - 1, v);
        const right = this.basisFunctionV(j + 1, m - 1, v);
        const denominator1 = this.knotsV[j + m] - this.knotsV[j];
        const denominator2 = this.knotsV[j + m + 1] - this.knotsV[j + 1];
        const term1 = denominator1 !== 0 ? (v - this.knotsV[j]) / denominator1 * left : 0;
        const term2 = denominator2 !== 0 ? (this.knotsV[j + m + 1] - v) / denominator2 * right : 0;
        return term1 + term2;
    }
    area(stepSizeU = 0.01, stepSizeV = 0.01) {
        let totalArea = 0;
        for (let u = 0; u < 1; u += stepSizeU) {
            for (let v = 0; v < 1; v += stepSizeV) {
                const p0 = this.evaluate(u, v);
                const p1 = this.evaluate(u + stepSizeU, v);
                const p2 = this.evaluate(u, v + stepSizeV);
                totalArea += this.triangleArea(p0, p1, p2);
            }
        }
        return totalArea;
    }
    trim(trimCurve) {
        const trimmedControlPoints = [];
        const trimmedKnotsU = [...this.knotsU];
        const trimmedKnotsV = [...this.knotsV];
        const controlPointSet = new Set();
        const trimPoints = [];
        const stepSizeU = (trimCurve.knots.length - 1) / 100;
        for (let u = trimCurve.knots[trimCurve.degree]; u <= trimCurve.knots[trimCurve.knots.length - trimCurve.degree - 1]; u += stepSizeU) {
            const trimPoint = trimCurve.evaluate(u);
            trimPoints.push(trimPoint);
        }
        for (const controlPointRow of this.controlPoints) {
            for (const controlPoint of controlPointRow) {
                if (this.isControlPointInsideTrimCurve(controlPoint, trimPoints)) {
                    controlPointSet.add(controlPoint);
                }
            }
        }
        for (const controlPointRow of this.controlPoints) {
            const row = controlPointRow.filter(cp => controlPointSet.has(cp));
            if (row.length > 0) {
                trimmedControlPoints.push(row);
            }
        }
        this.controlPoints = trimmedControlPoints;
        this.knotsU = trimmedKnotsU;
        this.knotsV = trimmedKnotsV;
    }
    intersection(otherSurface) {
        const intersections = [];
        const stepSizeU = 0.01;
        const stepSizeV = 0.01;
        for (let u = 0; u <= 1; u += stepSizeU) {
            for (let v = 0; v <= 1; v += stepSizeV) {
                const pointA = this.evaluate(u, v);
                for (let sU = 0; sU <= 1; sU += stepSizeU) {
                    for (let sV = 0; sV <= 1; sV += stepSizeV) {
                        const pointB = otherSurface.evaluate(sU, sV);
                        if (VectorUtils.isClose(pointA, pointB)) {
                            intersections.push({
                                point: pointA,
                                params: {
                                    u,
                                    v,
                                    sU,
                                    sV
                                }
                            });
                        }
                    }
                }
            }
        }
        return intersections;
    }
    fitToPoints(points) {
        const degreeU = this.degreeU;
        const degreeV = this.degreeV;
        const newControlPoints = Array.from({ length: points.length }, () => Array.from({ length: points[0].length }, () => new Array(points[0][0].length).fill(0)));
        const stepSizeU = 1 / (points.length - 1);
        const stepSizeV = 1 / (points[0].length - 1);
        for (let i = 0; i < points.length; i++) {
            for (let j = 0; j < points[0].length; j++) {
                const u = i * stepSizeU;
                const v = j * stepSizeV;
                const basisU = this.basisFunctionU(this.findSpanU(u), degreeU, u);
                const basisV = this.basisFunctionV(this.findSpanV(v), degreeV, v);
                for (let k = 0; k <= degreeU; k++) {
                    for (let l = 0; l <= degreeV; l++) {
                        newControlPoints[k][l] = newControlPoints[k][l].map((val, idx) => val + basisU[k] * basisV[l] * points[i][j][idx]);
                    }
                }
            }
        }
        this.controlPoints = newControlPoints;
    }
    continuityAnalysis() {
        const continuityLevels = {
            C0: [],
            C1: [],
            C2: []
        };
        const controlPointsU = this.controlPoints.length;
        const controlPointsV = this.controlPoints[0].length;
        // C0 Continuity: Check if the surfaces share boundary points
        for (let i = 0; i < controlPointsU - 1; i++) {
            for (let j = 0; j < controlPointsV - 1; j++) {
                const pointA = this.controlPoints[i][j];
                const pointB = this.controlPoints[i + 1][j];
                const pointC = this.controlPoints[i][j + 1];
                const pointD = this.controlPoints[i + 1][j + 1];
                if (this.isClose(pointA, pointB)) {
                    continuityLevels.C0.push({
                        points: [
                            pointA,
                            pointB
                        ],
                        type: 'C0'
                    });
                }
                if (this.isClose(pointA, pointC)) {
                    continuityLevels.C0.push({
                        points: [
                            pointA,
                            pointC
                        ],
                        type: 'C0'
                    });
                }
                if (this.isClose(pointB, pointD)) {
                    continuityLevels.C0.push({
                        points: [
                            pointB,
                            pointD
                        ],
                        type: 'C0'
                    });
                }
                if (this.isClose(pointC, pointD)) {
                    continuityLevels.C0.push({
                        points: [
                            pointC,
                            pointD
                        ],
                        type: 'C0'
                    });
                }
            }
        }
        // C1 Continuity: Check if the tangent vectors at boundary points are the same
        for (let i = 0; i < controlPointsU - 1; i++) {
            for (let j = 0; j < controlPointsV - 1; j++) {
                const pointA = this.controlPoints[i][j];
                const pointB = this.controlPoints[i + 1][j];
                const pointC = this.controlPoints[i][j + 1];
                const pointD = this.controlPoints[i + 1][j + 1];
                const tangentAB = this.calculateTangent(pointA, pointB);
                const tangentAC = this.calculateTangent(pointA, pointC);
                const tangentBD = this.calculateTangent(pointB, pointD);
                const tangentCD = this.calculateTangent(pointC, pointD);
                if (this.isClose(tangentAB, tangentAC)) {
                    continuityLevels.C1.push({
                        tangents: [
                            tangentAB,
                            tangentAC
                        ],
                        type: 'C1'
                    });
                }
                if (this.isClose(tangentAB, tangentBD)) {
                    continuityLevels.C1.push({
                        tangents: [
                            tangentAB,
                            tangentBD
                        ],
                        type: 'C1'
                    });
                }
                if (this.isClose(tangentAC, tangentCD)) {
                    continuityLevels.C1.push({
                        tangents: [
                            tangentAC,
                            tangentCD
                        ],
                        type: 'C1'
                    });
                }
            }
        }
        // C2 Continuity: Compare curvature vectors to check for curvature continuity
        for (let i = 0; i < controlPointsU - 2; i++) {
            for (let j = 0; j < controlPointsV - 2; j++) {
                const pointA = this.controlPoints[i][j];
                const pointB = this.controlPoints[i + 1][j];
                const pointC = this.controlPoints[i][j + 1];
                const pointD = this.controlPoints[i + 1][j + 1];
                const curvatureA = this.calculateCurvature(pointA, pointB, pointC);
                const curvatureB = this.calculateCurvature(pointB, pointD, pointC);
                if (this.isClose(curvatureA, curvatureB)) {
                    continuityLevels.C2.push({
                        curvatures: [
                            curvatureA,
                            curvatureB
                        ],
                        type: 'C2'
                    });
                }
            }
        }
        return continuityLevels;
    }
    findSpanU(u) {
        const n = this.knotsU.length - this.degreeU - 1;
        if (u >= this.knotsU[n + 1])
            return n;
        if (u <= this.knotsU[this.degreeU])
            return this.degreeU;
        let low = this.degreeU;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (u < this.knotsU[mid] || u >= this.knotsU[mid + 1]) {
            if (u < this.knotsU[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    findSpanV(v) {
        const n = this.knotsV.length - this.degreeV - 1;
        if (v >= this.knotsV[n + 1])
            return n;
        if (v <= this.knotsV[this.degreeV])
            return this.degreeV;
        let low = this.degreeV;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (v < this.knotsV[mid] || v >= this.knotsV[mid + 1]) {
            if (v < this.knotsV[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    static fromControlPoints(degreeU, degreeV, controlPoints) {
        const knotsU = Array(degreeU + 1).fill(0).concat(Array(controlPoints.length).fill(1));
        const knotsV = Array(degreeV + 1).fill(0).concat(Array(controlPoints[0].length).fill(1));
        return new NURBS_Surface(degreeU, degreeV, controlPoints, knotsU, knotsV);
    }
    evaluateSurfaceDerivative(u, v, direction) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const derivativePoint = new Array(this.controlPoints[0][0].length).fill(0);
        if (direction === 'u') {
            for (let i = 0; i <= this.degreeU; i++) {
                for (let j = 0; j <= this.degreeV; j++) {
                    const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                    const basisDerivativeU = this.basisFunctionU(spanU, this.degreeU - 1, u);
                    for (let k = 0; k < derivativePoint.length; k++) {
                        derivativePoint[k] += basisDerivativeU * basisV[j] * controlPoint[k];
                    }
                }
            }
        } else if (direction === 'v') {
            for (let i = 0; i <= this.degreeU; i++) {
                for (let j = 0; j <= this.degreeV; j++) {
                    const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                    const basisDerivativeV = this.basisFunctionV(spanV, this.degreeV - 1, v);
                    for (let k = 0; k < derivativePoint.length; k++) {
                        derivativePoint[k] += basisU[i] * basisDerivativeV * controlPoint[k];
                    }
                }
            }
        } else {
            throw new Error('Invalid direction specified. Use \'u\' or \'v\'.');
        }
        return derivativePoint;
    }
    parameterize() {
        const controlPointsU = this.controlPoints.length;
        const controlPointsV = this.controlPoints[0].length;
        const newKnotsU = Array(this.degreeU + 1).fill(0);
        for (let i = 0; i <= controlPointsU - this.degreeU; i++) {
            newKnotsU.push(i / (controlPointsU - this.degreeU));
        }
        newKnotsU.push(...Array(this.degreeU + 1).fill(1));
        const newKnotsV = Array(this.degreeV + 1).fill(0);
        for (let j = 0; j <= controlPointsV - this.degreeV; j++) {
            newKnotsV.push(j / (controlPointsV - this.degreeV));
        }
        newKnotsV.push(...Array(this.degreeV + 1).fill(1));
        this.knotsU = newKnotsU;
        this.knotsV = newKnotsV;
    }
    projectToSurface(point) {
        const stepSizeU = 0.01;
        const stepSizeV = 0.01;
        let minDistance = Infinity;
        let closestPoint = null;
        for (let u = 0; u <= 1; u += stepSizeU) {
            for (let v = 0; v <= 1; v += stepSizeV) {
                const surfacePoint = this.evaluate(u, v);
                const distance = this.distanceBetweenPoints(surfacePoint, point);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestPoint = surfacePoint;
                }
            }
        }
        return closestPoint;
    }
    getBoundaryCurves() {
        const boundaryCurves = [];
        const controlPointRows = this.controlPoints.length;
        const controlPointCols = this.controlPoints[0].length;
        // Top boundary curve
        const topCurve = [];
        for (let j = 0; j < controlPointCols; j++) {
            topCurve.push(this.controlPoints[0][j]);
        }
        boundaryCurves.push(new NURBS_Curve(this.degreeU, topCurve, this.knotsU));
        // Bottom boundary curve
        const bottomCurve = [];
        for (let j = 0; j < controlPointCols; j++) {
            bottomCurve.push(this.controlPoints[controlPointRows - 1][j]);
        }
        boundaryCurves.push(new NURBS_Curve(this.degreeU, bottomCurve, this.knotsU));
        // Left boundary curve
        const leftCurve = [];
        for (let i = 0; i < controlPointRows; i++) {
            leftCurve.push(this.controlPoints[i][0]);
        }
        boundaryCurves.push(new NURBS_Curve(this.degreeV, leftCurve, this.knotsV));
        // Right boundary curve
        const rightCurve = [];
        for (let i = 0; i < controlPointRows; i++) {
            rightCurve.push(this.controlPoints[i][controlPointCols - 1]);
        }
        boundaryCurves.push(new NURBS_Curve(this.degreeV, rightCurve, this.knotsV));
        return boundaryCurves;
    }
    getSurfaceInfo() {
        return {
            degreeU: this.degreeU,
            degreeV: this.degreeV,
            controlPointsCount: this.controlPoints.length
        };
    }
    trimCurve(curve) {
        const trimmedControlPoints = [];
        const trimmedKnotsU = [...this.knotsU];
        const trimmedKnotsV = [...this.knotsV];
        const controlPointSet = new Set();
        const trimPoints = [];
        const stepSizeU = (curve.knots.length - 1) / 100;
        for (let u = curve.knots[curve.degree]; u <= curve.knots[curve.knots.length - curve.degree - 1]; u += stepSizeU) {
            const trimPoint = curve.evaluate(u);
            trimPoints.push(trimPoint);
        }
        for (const controlPointRow of this.controlPoints) {
            for (const controlPoint of controlPointRow) {
                if (this.isControlPointInsideTrimCurve(controlPoint, trimPoints)) {
                    controlPointSet.add(controlPoint);
                }
            }
        }
        for (const controlPointRow of this.controlPoints) {
            const row = controlPointRow.filter(cp => controlPointSet.has(cp));
            if (row.length > 0) {
                trimmedControlPoints.push(row);
            }
        }
        this.controlPoints = trimmedControlPoints;
        this.knotsU = trimmedKnotsU;
        this.knotsV = trimmedKnotsV;
    }
    evaluateNormalAtPoint(u, v) {
        const normal = this.crossProduct(this.evaluateSurfaceDerivative(u, v, 'u'), this.evaluateSurfaceDerivative(u, v, 'v'));
        return VectorUtils.normalizeVector(normal);
    }
    getSurfaceParameterization() {
        const knotsU = this.knotsU;
        const knotsV = this.knotsV;
        const controlPointsU = this.controlPoints.length;
        const controlPointsV = this.controlPoints[0].length;
        const parameterization = {
            u: [],
            v: []
        };
        // Calculate parameterization for U
        for (let i = 0; i <= controlPointsU - this.degreeU; i++) {
            parameterization.u.push(i / (controlPointsU - this.degreeU));
        }
        parameterization.u.unshift(0);
        parameterization.u.push(1);
        // Calculate parameterization for V
        for (let j = 0; j <= controlPointsV - this.degreeV; j++) {
            parameterization.v.push(j / (controlPointsV - this.degreeV));
        }
        parameterization.v.unshift(0);
        parameterization.v.push(1);
        return parameterization;
    }
    surfaceArea(stepSizeU, stepSizeV) {
        let totalArea = 0;
        const uSteps = Math.ceil(1 / stepSizeU);
        const vSteps = Math.ceil(1 / stepSizeV);
        for (let i = 0; i < uSteps; i++) {
            for (let j = 0; j < vSteps; j++) {
                const u = i * stepSizeU;
                const v = j * stepSizeV;
                const p0 = this.evaluate(u, v);
                const p1 = this.evaluate(u + stepSizeU, v);
                const p2 = this.evaluate(u, v + stepSizeV);
                totalArea += this.triangleArea(p0, p1, p2);
            }
        }
        return totalArea;
    }
    isClosed() {
        const firstRow = this.controlPoints[0];
        const lastRow = this.controlPoints[this.controlPoints.length - 1];
        return firstRow[0] === lastRow[0] && firstRow[firstRow.length - 1] === lastRow[lastRow.length - 1];
    }
    static fromBoundaryCurves(curves) {
        const degreeU = 1;
        // Default degree for U
        const degreeV = 1;
        // Default degree for V
        const controlPoints = [];
        // Extract control points from curves
        curves.forEach(curve => {
            const points = curve.getPoints(1 / 100);
            // Sample points from the curve
            controlPoints.push(points);
        });
        // Prepare knots for U and V
        const knotsU = Array(degreeU + 1).fill(0).concat(Array(controlPoints.length).fill(1));
        const knotsV = Array(degreeV + 1).fill(0).concat(Array(controlPoints[0].length).fill(1));
        return new NURBS_Surface(degreeU, degreeV, controlPoints, knotsU, knotsV);
    }
    parameterizeSurface() {
        const controlPointsU = this.controlPoints.length;
        const controlPointsV = this.controlPoints[0].length;
        const newKnotsU = [];
        const newKnotsV = [];
        for (let i = 0; i <= this.degreeU; i++) {
            newKnotsU.push(0);
        }
        for (let i = 0; i <= controlPointsU - this.degreeU; i++) {
            newKnotsU.push(i / (controlPointsU - this.degreeU));
        }
        for (let i = 0; i <= this.degreeU; i++) {
            newKnotsU.push(1);
        }
        for (let j = 0; j <= this.degreeV; j++) {
            newKnotsV.push(0);
        }
        for (let j = 0; j <= controlPointsV - this.degreeV; j++) {
            newKnotsV.push(j / (controlPointsV - this.degreeV));
        }
        for (let j = 0; j <= this.degreeV; j++) {
            newKnotsV.push(1);
        }
        this.knotsU = newKnotsU;
        this.knotsV = newKnotsV;
    }
    crossProduct(v1, v2) {
        return [
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        ];
    }
    triangleArea(p0, p1, p2) {
        const a = this.distanceBetweenPoints(p0, p1);
        const b = this.distanceBetweenPoints(p1, p2);
        const c = this.distanceBetweenPoints(p2, p0);
        const s = (a + b + c) / 2;
        return Math.sqrt(s * (s - a) * (s - b) * (s - c));
    }
    distanceBetweenPoints(p1, p2) {
        return Math.sqrt(p1.reduce((sum, val, idx) => sum + (val - p2[idx]) ** 2, 0));
    }
    isControlPointInsideTrimCurve(controlPoint, trimPoints) {
        const [x, y] = controlPoint;
        const minX = Math.min(...trimPoints.map(p => p[0]));
        const maxX = Math.max(...trimPoints.map(p => p[0]));
        const minY = Math.min(...trimPoints.map(p => p[1]));
        const maxY = Math.max(...trimPoints.map(p => p[1]));
        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }
    isClose(point1, point2, tolerance = 0.000001) {
        return point1.every((val, index) => Math.abs(val - point2[index]) < tolerance);
    }
    calculateTangent(point1, point2) {
        return point2.map((val, idx) => val - point1[idx]);
    }
    evaluateSurfaceNormal(u, v) {
        const surfacePoint = this.evaluate(u, v);
        const normal = this.normal(u, v);
        return {
            surfacePoint,
            normal
        };
    }
    calculateCurvature(pointA, pointB, pointC) {
        const tangentAB = this.calculateTangent(pointA, pointB);
        const tangentAC = this.calculateTangent(pointA, pointC);
        const curvature = this.crossProduct(tangentAB, tangentAC);
        return curvature;
    }
    computeFirstDerivatives() {
        const derivatives = [];
        const stepSize = 0.01;
        for (let u = 0; u <= 1; u += stepSize) {
            for (let v = 0; v <= 1; v += stepSize) {
                const derivativeU = this.evaluateSurfaceDerivative(u, v, 'u');
                const derivativeV = this.evaluateSurfaceDerivative(u, v, 'v');
                derivatives.push({
                    u: derivativeU,
                    v: derivativeV
                });
            }
        }
        return derivatives;
    }
    validateKnotVectors() {
        const nU = this.controlPoints.length;
        const nV = this.controlPoints[0].length;
        const knotsMultiplicityU = this.knotsU.reduce((acc, knot) => {
            acc[knot] = (acc[knot] || 0) + 1;
            return acc;
        }, {});
        const knotsMultiplicityV = this.knotsV.reduce((acc, knot) => {
            acc[knot] = (acc[knot] || 0) + 1;
            return acc;
        }, {});
        const minKnotU = this.knotsU[0];
        const maxKnotU = this.knotsU[this.knotsU.length - 1];
        const minKnotV = this.knotsV[0];
        const maxKnotV = this.knotsV[this.knotsV.length - 1];
        if (this.knotsU.length < this.degreeU + nU + 1) {
            throw new Error('Knot vector U must have at least degreeU + number of control points + 1 elements.');
        }
        if (minKnotU !== 0 || maxKnotU !== 1) {
            throw new Error('Knot vector U must start with 0 and end with 1.');
        }
        for (const knot in knotsMultiplicityU) {
            if (knotsMultiplicityU[knot] > this.degreeU + 1) {
                throw new Error(`Knot ${ knot } in U has a multiplicity greater than degreeU + 1.`);
            }
        }
        if (this.knotsV.length < this.degreeV + nV + 1) {
            throw new Error('Knot vector V must have at least degreeV + number of control points + 1 elements.');
        }
        if (minKnotV !== 0 || maxKnotV !== 1) {
            throw new Error('Knot vector V must start with 0 and end with 1.');
        }
        for (const knot in knotsMultiplicityV) {
            if (knotsMultiplicityV[knot] > this.degreeV + 1) {
                throw new Error(`Knot ${ knot } in V has a multiplicity greater than degreeV + 1.`);
            }
        }
    }
    reorderKnotVectors(newKnotsU, newKnotsV) {
        this.knotsU = newKnotsU;
        this.knotsV = newKnotsV;
        this.validateKnotVectors();
    }
    normalizeParameterization() {
        const n = this.controlPoints.length - 1;
        const newKnots = [];
        const newControlPoints = [...this.controlPoints];
        const totalLength = this.length();
        if (totalLength === 0) {
            throw new Error('Total length of the curve is zero. Cannot normalize.');
        }
        const accumulatedLength = new Array(n + 1).fill(0);
        for (let i = 0; i < n; i++) {
            const p0 = this.evaluate(i / n);
            const p1 = this.evaluate((i + 1) / n);
            accumulatedLength[i + 1] = accumulatedLength[i] + this.distanceBetweenPoints(p0, p1);
        }
        for (let i = 0; i <= n; i++) {
            const normalizedParameter = accumulatedLength[i] / totalLength;
            newKnots.push(i === 0 || i === n ? normalizedParameter : normalizedParameter);
        }
        this.knots = newKnots;
        this.controlPoints = newControlPoints;
    }
    generateControlPoints() {
        const controlPoints = [];
        const numSegments = 36;
        const numStacks = 18;
        for (let i = 0; i <= numStacks; i++) {
            const stackAngle = Math.PI / 2 - i * Math.PI / numStacks;
            const xy = this.radius * Math.cos(stackAngle);
            const z = this.radius * Math.sin(stackAngle);
            const row = [];
            for (let j = 0; j <= numSegments; j++) {
                const sectorAngle = j * 2 * Math.PI / numSegments;
                const x = this.center[0] + xy * Math.cos(sectorAngle);
                const y = this.center[1] + xy * Math.sin(sectorAngle);
                row.push([
                    x,
                    y,
                    z + this.center[2]
                ]);
            }
            controlPoints.push(row);
        }
        this.controlPoints = controlPoints;
    }
    toBufferGeometry(stepSizeU = 0.1, stepSizeV = 0.1) {
        const points = this.getPoints(stepSizeU, stepSizeV);
        const geometry = new THREE.BufferGeometry();
        const vertices = new Float32Array(points.flat());
        const indices = [];
        const width = Math.ceil(1 / stepSizeU) + 1;
        for (let i = 0; i < this.controlPoints.length - 1; i++) {
            for (let j = 0; j < width - 1; j++) {
                const a = i * width + j;
                const b = a + width;
                indices.push(a, b, a + 1);
                indices.push(b, b + 1, a + 1);
            }
        }
        geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
        geometry.setIndex(indices);
        return geometry;
    }
    intersect(otherSurface) {
        const intersections = [];
        const stepSizeU = 0.01;
        const stepSizeV = 0.01;
        for (let u = 0; u <= 1; u += stepSizeU) {
            for (let v = 0; v <= 1; v += stepSizeV) {
                const pointA = this.evaluate(u, v);
                for (let sU = 0; sU <= 1; sU += stepSizeU) {
                    for (let sV = 0; sV <= 1; sV += stepSizeV) {
                        const pointB = otherSurface.evaluate(sU, sV);
                        if (VectorUtils.isClose(pointA, pointB)) {
                            intersections.push({
                                point: pointA,
                                params: {
                                    u,
                                    v,
                                    sU,
                                    sV
                                }
                            });
                        }
                    }
                }
            }
        }
        return intersections;
    }
    intersectsWith(otherSurface) {
        const intersections = [];
        const stepSizeU = 0.01;
        const stepSizeV = 0.01;
        for (let u = 0; u <= 1; u += stepSizeU) {
            for (let v = 0; v <= 1; v += stepSizeV) {
                const pointA = this.evaluate(u, v);
                for (let sU = 0; sU <= 1; sU += stepSizeU) {
                    for (let sV = 0; sV <= 1; sV += stepSizeV) {
                        const pointB = otherSurface.evaluate(sU, sV);
                        if (VectorUtils.isClose(pointA, pointB)) {
                            intersections.push({
                                point: pointA,
                                paramsA: {
                                    u,
                                    v
                                },
                                paramsB: {
                                    sU,
                                    sV
                                }
                            });
                        }
                    }
                }
            }
        }
        return intersections;
    }
    closestIntersection(otherSurface) {
        const intersections = this.intersectsWith(otherSurface);
        if (intersections.length === 0)
            return null;
        let minDistance = Infinity;
        let closestIntersection = null;
        intersections.forEach(({point, paramsA, paramsB}) => {
            const distance = VectorUtils.distanceBetweenPoints(point, this.evaluate(paramsA.u, paramsA.v));
            if (distance < minDistance) {
                minDistance = distance;
                closestIntersection = {
                    point,
                    paramsA,
                    paramsB
                };
            }
        });
        return closestIntersection;
    }
    parameterizeKnots(n, degree) {
        const knots = [];
        for (let i = 0; i <= degree; i++) {
            knots.push(0);
        }
        for (let i = 0; i < n - degree; i++) {
            knots.push(i / (n - degree));
        }
        for (let i = 0; i <= degree; i++) {
            knots.push(1);
        }
        return knots;
    }
    isPointInsideTrimCurve(point, trimPoints) {
        const [x, y] = point;
        const minX = Math.min(...trimPoints.map(p => p[0]));
        const maxX = Math.max(...trimPoints.map(p => p[0]));
        const minY = Math.min(...trimPoints.map(p => p[1]));
        const maxY = Math.max(...trimPoints.map(p => p[1]));
        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }
    evaluateSurfaceAtMultipleParameters(parameterPairs) {
        return parameterPairs.map(([u, v]) => this.evaluate(u, v));
    }
}
class Line extends NURBS_Curve {
    constructor(startPoint, endPoint) {
        const controlPoints = [
            [
                ...startPoint,
                0
            ],
            // Ensure 3D
            [
                ...endPoint,
                0
            ]
        ];
        // Ensure 3D
        const degree = 1;
        const knots = [
            0,
            0,
            1,
            1
        ];
        super(degree, controlPoints, knots);
    }
    getMidPoint() {
        return this.controlPoints[0].map((val, index) => (val + this.controlPoints[1][index]) / 2);
    }
}
class Circle extends NURBS_Curve {
    constructor(center, radius) {
        super(2, [], []);
        this.center = center;
        this.radius = radius;
    }
    getCircumference() {
        return 2 * Math.PI * this.radius;
    }
    generateControlPoints() {
        const numPoints = 100;
        for (let i = 0; i < numPoints; i++) {
            const angle = i / numPoints * 2 * Math.PI;
            const x = this.center[0] + this.radius * Math.cos(angle);
            const y = this.center[1] + this.radius * Math.sin(angle);
            this.controlPoints.push([
                x,
                y,
                this.center[2]
            ]);
        }
    }
    getArea() {
        return Math.PI * this.radius ** 2;
    }
    getPointAtAngle(angle) {
        const x = this.center[0] + this.radius * Math.cos(angle);
        const y = this.center[1] + this.radius * Math.sin(angle);
        return [
            x,
            y,
            this.center[2]
        ];
    }
    getTangentAtAngle(angle) {
        const tangent = [
            -this.radius * Math.sin(angle),
            this.radius * Math.cos(angle),
            0
        ];
        const length = Math.sqrt(tangent[0] ** 2 + tangent[1] ** 2);
        return tangent.map(val => val / length);
    }
}
// Using radius value from control points
class Arc extends NURBS_Curve {
    constructor(center, radius, startAngle, endAngle) {
        super(2, [], []);
        this.center = center;
        this.radius = radius;
        this.startAngle = startAngle;
        this.endAngle = endAngle;
    }
    getLength() {
        const angleSweep = this.endAngle - this.startAngle;
        return Math.abs(angleSweep) * this.radius;
    }
    getCenter() {
        const cx = this.center[0];
        const cy = this.center[1];
        return [
            cx,
            cy
        ];
    }
    generateControlPoints() {
        const numPoints = 100;
        for (let i = 0; i <= numPoints; i++) {
            const angle = this.startAngle + i / numPoints * (this.endAngle - this.startAngle);
            const x = this.center[0] + this.radius * Math.cos(angle);
            const y = this.center[1] + this.radius * Math.sin(angle);
            this.controlPoints.push([
                x,
                y,
                this.center[2]
            ]);
        }
    }
    getPointAtAngle(angle) {
        const x = this.center[0] + this.radius * Math.cos(angle);
        const y = this.center[1] + this.radius * Math.sin(angle);
        return [
            x,
            y,
            this.center[2]
        ];
    }
    getTangentAtAngle(angle) {
        const tangent = [
            -this.radius * Math.sin(angle),
            this.radius * Math.cos(angle),
            0
        ];
        const length = Math.sqrt(tangent[0] ** 2 + tangent[1] ** 2);
        return tangent.map(val => val / length);
    }
}
class Ellipse extends NURBS_Curve {
    constructor(center, radiusX, radiusY) {
        super(2, [], []);
        this.center = center;
        this.radiusX = radiusX;
        this.radiusY = radiusY;
    }
    getCircumference() {
        return 2 * Math.PI * Math.sqrt((this.radiusX ** 2 + this.radiusY ** 2) / 2);
    }
    getArea() {
        return Math.PI * this.radiusX * this.radiusY;
    }
}
class BsplinesCurve extends NURBS_Curve {
    constructor(controlPoints, degree) {
        if (controlPoints.length < degree + 1) {
            throw new Error(`BSpline curve requires at least ${ degree + 1 } control points.`);
        }
        const n = controlPoints.length;
        const knots = Array(degree + 1).fill(0).concat(Array.from({ length: n - degree - 1 }, (_, i) => i / (n - degree - 1)), Array(degree + 1).fill(1));
        super(degree, controlPoints, knots);
    }
    getDegree() {
        return this.degree;
    }
    getControlPoints() {
        return this.controlPoints.slice();
    }
}
class EllipticalArc extends NURBS_Curve {
    constructor(center, radiusX, radiusY, startAngle, endAngle) {
        super(2, [], []);
        this.center = center;
        this.radiusX = radiusX;
        this.radiusY = radiusY;
        this.startAngle = startAngle;
        this.endAngle = endAngle;
        this.generateControlPoints();
    }
    getLength() {
        const theta1 = this.startAngle;
        const theta2 = this.endAngle;
        const arcLength = (theta2 - theta1) * Math.sqrt(this.radiusX ** 2 + this.radiusY ** 2) / 2;
        return Math.abs(arcLength);
    }
    generateControlPoints() {
        const numPoints = 100;
        for (let i = 0; i <= numPoints; i++) {
            const angle = this.startAngle + i / numPoints * (this.endAngle - this.startAngle);
            const x = this.center[0] + this.radiusX * Math.cos(angle);
            const y = this.center[1] + this.radiusY * Math.sin(angle);
            this.controlPoints.push([
                x,
                y,
                this.center[2]
            ]);
        }
    }
    getCenter() {
        return this.center;
    }
    getPointAtAngle(angle) {
        const x = this.center[0] + this.radiusX * Math.cos(angle);
        const y = this.center[1] + this.radiusY * Math.sin(angle);
        return [
            x,
            y,
            this.center[2]
        ];
    }
}
class QuadraticBezierCurve extends NURBS_Curve {
    constructor(controlPoints) {
        super(2, controlPoints, []);
    }
    getQuadraticBezierPoint(t) {
        if (t < 0 || t > 1) {
            throw new Error('Parameter t must be between 0 and 1.');
        }
        const [P0, P1, P2] = this.controlPoints;
        const x = (1 - t) ** 2 * P0[0] + 2 * (1 - t) * t * P1[0] + t ** 2 * P2[0];
        const y = (1 - t) ** 2 * P0[1] + 2 * (1 - t) * t * P1[1] + t ** 2 * P2[1];
        return [
            x,
            y,
            P0[2]
        ];
    }
    // Assuming Z coordinate remains constant
    getTangent(t) {
        if (t < 0 || t > 1) {
            throw new Error('Parameter t must be between 0 and 1.');
        }
        const [P0, P1, P2] = this.controlPoints;
        const tangentX = 2 * (1 - t) * (P1[0] - P0[0]) + 2 * t * (P2[0] - P1[0]);
        const tangentY = 2 * (1 - t) * (P1[1] - P0[1]) + 2 * t * (P2[1] - P1[1]);
        const length = Math.sqrt(tangentX ** 2 + tangentY ** 2);
        if (length === 0) {
            throw new Error('The tangent at the parameter t results in a zero vector. Cannot compute tangent.');
        }
        return [
            tangentX / length,
            tangentY / length,
            0
        ];
    }
}
class CubicBezierCurve extends NURBS_Curve {
    constructor(controlPoints) {
        super(3, controlPoints, []);
    }
    getCubicBezierPoint(t) {
        if (t < 0 || t > 1) {
            throw new Error('Parameter t must be between 0 and 1.');
        }
        const [P0, P1, P2, P3] = this.controlPoints;
        const u = 1 - t;
        const x = u ** 3 * P0[0] + 3 * u ** 2 * t * P1[0] + 3 * u * t ** 2 * P2[0] + t ** 3 * P3[0];
        const y = u ** 3 * P0[1] + 3 * u ** 2 * t * P1[1] + 3 * u * t ** 2 * P2[1] + t ** 3 * P3[1];
        return [
            x,
            y,
            P0[2]
        ];
    }
    getTangent(t) {
        if (t < 0 || t > 1) {
            throw new Error('Parameter t must be between 0 and 1.');
        }
        const [P0, P1, P2, P3] = this.controlPoints;
        const u = 1 - t;
        const tangentX = -3 * u ** 2 * P0[0] + 3 * (3 * u ** 2 - 2 * u * t) * P1[0] + 3 * (2 * u * t - 3 * t ** 2) * P2[0] + 3 * t ** 2 * P3[0];
        const tangentY = -3 * u ** 2 * P0[1] + 3 * (3 * u ** 2 - 2 * u * t) * P1[1] + 3 * (2 * u * t - 3 * t ** 2) * P2[1] + 3 * t ** 2 * P3[1];
        const length = Math.sqrt(tangentX ** 2 + tangentY ** 2);
        if (length === 0) {
            throw new Error('The tangent at the parameter t results in a zero vector. Cannot compute tangent.');
        }
        return [
            tangentX / length,
            tangentY / length,
            0
        ];
    }
}
class Parabola extends NURBS_Curve {
    constructor(vertex, focus) {
        super(2, [
            vertex,
            focus
        ], []);
    }
    getFocus() {
        return this.controlPoints[1];
    }
    getVertex() {
        return this.controlPoints[0];
    }
    getDirectrix() {
        const vertex = this.getVertex();
        const focus = this.getFocus();
        const directrixY = (vertex[1] + focus[1]) / 2;
        return directrixY;
    }
    isPointOnParabola(point) {
        const vertex = this.getVertex();
        const focus = this.getFocus();
        const directrixY = this.getDirectrix();
        const distanceToFocus = this.distanceBetweenPoints(point, focus);
        const distanceToDirectrix = Math.abs(point[1] - directrixY);
        return Math.abs(distanceToFocus - distanceToDirectrix) < 0.000001;
    }
}
class Hyperbola extends NURBS_Curve {
    constructor(focus1, focus2) {
        super(2, [
            focus1,
            focus2
        ], []);
    }
    getAsymptotes() {
        const [focus1, focus2] = this.controlPoints;
        const center = this.getCenter();
        const distance = this.distanceBetweenPoints(focus1, focus2) / 2;
        const slope1 = (focus1[1] - focus2[1]) / (focus1[0] - focus2[0]);
        const slope2 = -1 / slope1;
        const asymptote1 = {
            point: center,
            slope: slope1
        };
        const asymptote2 = {
            point: center,
            slope: slope2
        };
        return [
            asymptote1,
            asymptote2
        ];
    }
    getVertices() {
        const [focus1, focus2] = this.controlPoints;
        const center = this.getCenter();
        const distance = this.getDistanceBetweenFoci() / 2;
        return [
            [
                center[0] - distance,
                center[1],
                center[2]
            ],
            [
                center[0] + distance,
                center[1],
                center[2]
            ]
        ];
    }
    isPointOnHyperbola(point) {
        const vertices = this.getVertices();
        const distanceToFocus1 = this.distanceBetweenPoints(point, this.controlPoints[0]);
        const distanceToFocus2 = this.distanceBetweenPoints(point, this.controlPoints[1]);
        const distanceBetweenFoci = this.getDistanceBetweenFoci();
        return Math.abs(distanceToFocus1 - distanceToFocus2) === distanceBetweenFoci;
    }
    getDistanceBetweenFoci() {
        return this.distanceBetweenPoints(this.controlPoints[0], this.controlPoints[1]);
    }
    getEccentricity() {
        const distance = this.getDistanceBetweenFoci();
        const vertices = this.getVertices();
        const distanceVertices = this.distanceBetweenPoints(vertices[0], vertices[1]);
        return distance / distanceVertices;
    }
    getCenter() {
        const [focus1, focus2] = this.controlPoints;
        return [
            (focus1[0] + focus2[0]) / 2,
            (focus1[1] + focus2[1]) / 2,
            (focus1[2] + focus2[2]) / 2
        ];
    }
}
class Spiral extends NURBS_Curve {
    constructor(center, startAngle, growthFactor) {
        super(2, [], []);
        this.center = center;
        this.startAngle = startAngle;
        this.growthFactor = growthFactor;
        this.generateControlPoints();
    }
    // Add new control point
    getLength(stepSize = 0.01) {
        let totalLength = 0;
        let previousPoint = this.controlPoints[0];
        for (let i = 1; i < this.controlPoints.length; i++) {
            const currentPoint = this.controlPoints[i];
            totalLength += this.distanceBetweenPoints(previousPoint, currentPoint);
            previousPoint = currentPoint;
        }
        return totalLength;
    }
    generateControlPoints() {
        const numPoints = 100;
        for (let i = 0; i <= numPoints; i++) {
            const angle = this.startAngle + i * this.growthFactor;
            const radius = i / numPoints;
            // Adjust radius for spiral effect
            const x = this.center[0] + radius * Math.cos(angle);
            const y = this.center[1] + radius * Math.sin(angle);
            this.controlPoints.push([
                x,
                y,
                this.center[2]
            ]);
        }
    }
    getPointAtAngle(angle) {
        const radius = this.growthFactor * angle;
        // Calculate radius based on angle
        const x = this.center[0] + radius * Math.cos(angle);
        const y = this.center[1] + radius * Math.sin(angle);
        return [
            x,
            y,
            this.center[2]
        ];
    }
}
class PolynomialCurve extends NURBS_Curve {
    constructor(controlPoints, degree) {
        super(degree, controlPoints, []);
    }
    getValueAt(t) {
        return this.controlPoints.reduce((acc, pt, i) => acc + pt * Math.pow(t, i), 0);
    }
    getDerivativeAt(t) {
        return this.controlPoints.reduce((acc, pt, i) => acc + pt * i * Math.pow(t, i - 1), 0);
    }
}
class CatmullRomSpline extends NURBS_Curve {
    constructor(controlPoints) {
        super(2, controlPoints, []);
    }
    getPoint(t) {
        const p = this.controlPoints;
        const n = p.length;
        if (n < 4) {
            throw new Error('At least 4 control points are required.');
        }
        const t2 = t * t;
        const t3 = t2 * t;
        const p0 = p[Math.floor(t) % n];
        const p1 = p[(Math.floor(t) + 1) % n];
        const p2 = p[(Math.floor(t) + 2) % n];
        const p3 = p[(Math.floor(t) + 3) % n];
        return [
            0.5 * (2 * p1[0] + (-p0[0] + p2[0]) * t + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t2 + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t3),
            0.5 * (2 * p1[1] + (-p0[1] + p2[1]) * t + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t2 + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t3)
        ];
    }
    getTangent(t) {
        const p = this.controlPoints;
        const n = p.length;
        if (n < 4) {
            throw new Error('At least 4 control points are required.');
        }
        const t2 = t * t;
        const p0 = p[Math.floor(t) % n];
        const p1 = p[(Math.floor(t) + 1) % n];
        const p2 = p[(Math.floor(t) + 2) % n];
        const p3 = p[(Math.floor(t) + 3) % n];
        return [
            0.5 * (-p0[0] + p2[0] + (2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]) * t + (-p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]) * t2),
            0.5 * (-p0[1] + p2[1] + (2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]) * t + (-p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]) * t2)
        ];
    }
}
class HermiteCurve extends NURBS_Curve {
    constructor(pointA, tangentA, pointB, tangentB) {
        super(2, [
            pointA,
            pointB
        ], []);
        this.tangentA = tangentA;
        this.tangentB = tangentB;
    }
    getPoint(t) {
        if (t < 0 || t > 1) {
            throw new Error('Parameter t must be between 0 and 1.');
        }
        const [P0, P1] = this.controlPoints;
        const u = 1 - t;
        const point = [
            u ** 3 * P0[0] + 3 * u ** 2 * t * this.tangentA[0] + 3 * u * t ** 2 * this.tangentB[0] + t ** 3 * P1[0],
            u ** 3 * P0[1] + 3 * u ** 2 * t * this.tangentA[1] + 3 * u * t ** 2 * this.tangentB[1] + t ** 3 * P1[1],
            P0[2]
        ];
        return point;
    }
    getTangent(t) {
        if (t < 0 || t > 1) {
            throw new Error('Parameter t must be between 0 and 1.');
        }
        const [P0, P1] = this.controlPoints;
        const u = 1 - t;
        const tangent = [
            3 * u ** 2 * (P1[0] - P0[0]) + 6 * u * t * this.tangentA[0] + 3 * t ** 2 * this.tangentB[0],
            3 * u ** 2 * (P1[1] - P0[1]) + 6 * u * t * this.tangentA[1] + 3 * t ** 2 * this.tangentB[1],
            0
        ];
        const length = Math.sqrt(tangent[0] ** 2 + tangent[1] ** 2);
        if (length === 0) {
            throw new Error('The tangent at the parameter t results in a zero vector. Cannot compute tangent.');
        }
        return tangent.map(val => val / length);
    }
}
class Plane extends NURBS_Surface {
    constructor(point, normal) {
        super(1, 1, [], [], []);
        this.point = point;
        this.normal = this.normalize(normal);
        this.generateControlPoints();
    }
    generateControlPoints() {
        const uSteps = 5;
        const vSteps = 5;
        this.controlPoints = Array.from({ length: uSteps + 1 }, (_, i) => {
            return Array.from({ length: vSteps + 1 }, (_, j) => {
                return [
                    this.point[0] + i / uSteps * this.normal[0],
                    this.point[1] + j / vSteps * this.normal[1],
                    this.point[2]
                ];
            });
        });
    }
    // Keeping z constant for a simple plane
    evaluate(u, v) {
        const [spanU, spanV] = [
            this.findSpanU(u),
            this.findSpanV(v)
        ];
        const [basisU, basisV] = [
            this.basisFunctionU(spanU, this.degreeU, u),
            this.basisFunctionV(spanV, this.degreeV, v)
        ];
        const point = new Array(this.controlPoints[0][0].length).fill(0);
        for (let i = 0; i <= this.degreeU; i++) {
            for (let j = 0; j <= this.degreeV; j++) {
                const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                for (let k = 0; k < point.length; k++) {
                    point[k] += basisU[i] * basisV[j] * controlPoint[k];
                }
            }
        }
        return point;
    }
    normalize(vector) {
        const length = Math.sqrt(vector.reduce((sum, val) => sum + val ** 2, 0));
        if (length === 0)
            throw new Error('Cannot normalize a zero-length vector.');
        return vector.map(val => val / length);
    }
}
class Sphere extends NURBS_Surface {
    constructor(center, radius) {
        super(2, 2, [], [], []);
        this.center = center;
        this.radius = radius;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const controlPoints = [];
        const numSegments = 36;
        // Control points in radial direction
        const numStacks = 18;
        // Control points in vertical direction
        for (let i = 0; i <= numStacks; i++) {
            const stackAngle = Math.PI / 2 - i * Math.PI / numStacks;
            // from pi/2 to -pi/2
            const xy = this.radius * Math.cos(stackAngle);
            const z = this.radius * Math.sin(stackAngle);
            const row = [];
            for (let j = 0; j <= numSegments; j++) {
                const sectorAngle = j * 2 * Math.PI / numSegments;
                // from 0 to 2pi
                const x = this.center[0] + xy * Math.cos(sectorAngle);
                const y = this.center[1] + xy * Math.sin(sectorAngle);
                row.push([
                    x,
                    y,
                    this.center[2] + z
                ]);
            }
            controlPoints.push(row);
        }
        this.controlPoints = controlPoints;
        this.knotsU = Array(this.degreeU + 1).fill(0).concat(Array(numStacks - this.degreeU).fill(1));
        this.knotsV = Array(this.degreeV + 1).fill(0).concat(Array(numSegments - this.degreeV).fill(1));
    }
    evaluate(u, v) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const point = new Array(3).fill(0);
        for (let i = 0; i <= this.degreeU; i++) {
            for (let j = 0; j <= this.degreeV; j++) {
                const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                for (let k = 0; k < point.length; k++) {
                    point[k] += basisU[i] * basisV[j] * controlPoint[k];
                }
            }
        }
        return point;
    }
}
class Cylinder extends NURBS_Surface {
    constructor(axisPoint, axisDirection, radius, height) {
        super(2, 1, [], [], []);
        this.axisPoint = axisPoint;
        this.axisDirection = this.normalize(axisDirection);
        this.radius = radius;
        this.height = height;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const points = [];
        const numSegments = 36;
        const axisBase = this.axisPoint;
        for (let i = 0; i <= numSegments; i++) {
            const angle = i / numSegments * 2 * Math.PI;
            const x = axisBase[0] + this.radius * Math.cos(angle);
            const y = axisBase[1] + this.radius * Math.sin(angle);
            const zTop = axisBase[2] + this.height;
            // Top circle
            const zBottom = axisBase[2];
            // Bottom circle
            // Add control points for bottom and top circles
            points.push([
                x,
                y,
                zBottom
            ]);
            points.push([
                x,
                y,
                zTop
            ]);
        }
        this.controlPoints = points;
        this.knotsU = Array(2).fill(0).concat(Array(numSegments + 1).fill(1));
        this.knotsV = [
            0,
            0,
            1,
            1
        ];
    }
    evaluate(u, v) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const point = [
            0,
            0,
            0
        ];
        // Initialize point
        for (let i = 0; i <= this.degreeU; i++) {
            for (let j = 0; j <= this.degreeV; j++) {
                const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                for (let k = 0; k < point.length; k++) {
                    point[k] += basisU[i] * basisV[j] * controlPoint[k];
                }
            }
        }
        return point;
    }
    normalize(vector) {
        const length = Math.sqrt(vector.reduce((sum, val) => sum + val ** 2, 0));
        if (length === 0) {
            throw new Error('Cannot normalize a zero-length vector.');
        }
        return vector.map(val => val / length);
    }
}
class Cone extends NURBS_Surface {
    constructor(vertex, axisDirection, baseRadius, height) {
        super(2, 1, [], [], []);
        this.vertex = vertex;
        this.axisDirection = this.normalize(axisDirection);
        this.baseRadius = baseRadius;
        this.height = height;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const points = [];
        const numSegments = 36;
        for (let i = 0; i <= numSegments; i++) {
            const angle = i / numSegments * 2 * Math.PI;
            const x = this.baseRadius * Math.cos(angle);
            const y = this.baseRadius * Math.sin(angle);
            const z = this.vertex[2];
            // Bottom vertex
            points.push([
                x + this.vertex[0],
                y + this.vertex[1],
                z
            ]);
        }
        points.push([
            this.vertex[0],
            this.vertex[1],
            this.vertex[2] + this.height
        ]);
        // Apex of the cone
        this.controlPoints = points;
        this.knotsU = Array(2).fill(0).concat(Array(numSegments).fill(1));
        this.knotsV = [
            0,
            0,
            1,
            1
        ];
    }
    evaluate(u, v) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const point = [
            0,
            0,
            0
        ];
        for (let i = 0; i <= this.degreeU; i++) {
            for (let j = 0; j <= this.degreeV; j++) {
                const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                for (let k = 0; k < point.length; k++) {
                    point[k] += basisU[i] * basisV[j] * controlPoint[k];
                }
            }
        }
        return point;
    }
}
class Torus extends NURBS_Surface {
    constructor(center, majorRadius, minorRadius) {
        super(2, 2, [], [], []);
        this.center = center;
        this.majorRadius = majorRadius;
        this.minorRadius = minorRadius;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const controlPoints = [];
        const numSegments = 36;
        const numRadial = 18;
        for (let i = 0; i <= numRadial; i++) {
            const v = i / numRadial * 2 * Math.PI;
            // From 0 to 2π for circular motion
            const z = this.minorRadius * Math.sin(v);
            // Y-component
            const row = [];
            for (let j = 0; j <= numSegments; j++) {
                const u = j / numSegments * 2 * Math.PI;
                // From 0 to 2π for tubular motion
                const x = (this.majorRadius + z) * Math.cos(u);
                const y = (this.majorRadius + z) * Math.sin(u);
                row.push([
                    x + this.center[0],
                    y + this.center[1],
                    z + this.center[2]
                ]);
            }
            controlPoints.push(row);
        }
        this.controlPoints = controlPoints;
    }
    evaluate(u, v) {
        const spanU = this.findSpanU(u);
        const spanV = this.findSpanV(v);
        const basisU = this.basisFunctionU(spanU, this.degreeU, u);
        const basisV = this.basisFunctionV(spanV, this.degreeV, v);
        const point = new Array(3).fill(0);
        for (let i = 0; i <= this.degreeU; i++) {
            for (let j = 0; j <= this.degreeV; j++) {
                const controlPoint = this.controlPoints[spanU - this.degreeU + i][spanV - this.degreeV + j];
                for (let k = 0; k < point.length; k++) {
                    point[k] += basisU[i] * basisV[j] * controlPoint[k];
                }
            }
        }
        return point;
    }
}
class Box extends NURBS_Surface {
    constructor(minCorner, maxCorner) {
        super(1, 1, [], [], []);
        this.minCorner = minCorner;
        this.maxCorner = maxCorner;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const controlPoints = [];
        const corners = [
            this.minCorner,
            [
                this.maxCorner[0],
                this.minCorner[1],
                this.minCorner[2]
            ],
            this.maxCorner,
            [
                this.minCorner[0],
                this.maxCorner[1],
                this.maxCorner[2]
            ]
        ];
        for (const corner of corners) {
            controlPoints.push(corner);
        }
        this.controlPoints = [controlPoints];
    }
    evaluate(u, v) {
        const x = this.minCorner[0] + u * (this.maxCorner[0] - this.minCorner[0]);
        const y = this.minCorner[1] + v * (this.maxCorner[1] - this.minCorner[1]);
        const z = (this.minCorner[2] + this.maxCorner[2]) / 2;
        return [
            x,
            y,
            z
        ];
    }
}
class BezierSurface extends NURBS_Surface {
    constructor(controlPoints) {
        super(controlPoints.length - 1, controlPoints[0].length - 1, controlPoints, [], []);
        this.generateControlPoints();
    }
    // In this case, it's direct assignment from the constructor.
    evaluate(u, v) {
        const degreeU = this.degreeU;
        const degreeV = this.degreeV;
        const point = [];
        for (let i = 0; i <= degreeU; i++) {
            for (let j = 0; j <= degreeV; j++) {
                const controlPoint = this.controlPoints[i][j];
                const basisU = this.basisFunctionU(i, degreeU, u);
                const basisV = this.basisFunctionV(j, degreeV, v);
                point.forEach((val, index) => {
                    point[index] += basisU * basisV * controlPoint[index];
                });
            }
        }
        return point;
    }
    generateControlPoints() {
        const controlPoints = this.controlPoints;
        this.controlPoints = controlPoints;
    }
}
class RuledSurface extends NURBS_Surface {
    constructor(curve1, curve2) {
        super(1, 1, [], [], []);
        this.curve1 = curve1;
        this.curve2 = curve2;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const controlPoints = [];
        const segments = 20;
        // Number of segments for the ruled surface
        for (let i = 0; i <= segments; i++) {
            const u = i / segments;
            const point1 = this.curve1.evaluate(u);
            const point2 = this.curve2.evaluate(u);
            controlPoints.push([
                point1,
                point2
            ]);
        }
        this.controlPoints = controlPoints;
    }
    evaluate(u, v) {
        const point1 = this.curve1.evaluate(u);
        const point2 = this.curve2.evaluate(u);
        return point1.map((coord, index) => coord * (1 - v) + point2[index] * v);
    }
}
class LoftedSurface extends NURBS_Surface {
    constructor(curves) {
        super(curves[0].degree, curves.length - 1, [], [], []);
        this.curves = curves;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const controlPoints = [];
        const segments = 20;
        // Number of segments for the lofted surface
        for (let i = 0; i < this.curves.length; i++) {
            const curve = this.curves[i];
            const curvePoints = [];
            for (let j = 0; j <= segments; j++) {
                const u = j / segments;
                const point = curve.evaluate(u);
                curvePoints.push(point);
            }
            controlPoints.push(curvePoints);
        }
        this.controlPoints = controlPoints;
    }
    evaluate(u, v) {
        const idx = Math.floor(v * (this.curves.length - 1));
        const nextIdx = Math.min(idx + 1, this.curves.length - 1);
        const point1 = this.curves[idx].evaluate(u);
        const point2 = this.curves[nextIdx].evaluate(u);
        return point1.map((coord, index) => coord * (1 - (v - idx)) + point2[index] * (v - idx));
    }
}
class RevolutionSurface extends NURBS_Surface {
    constructor(profileCurve, axis) {
        super(1, 1, [], [], []);
        this.profileCurve = profileCurve;
        this.axis = axis;
        this.generateControlPoints();
    }
    generateControlPoints() {
        const controlPoints = [];
        const segments = 36;
        // Number of segments for the revolution
        for (let i = 0; i <= segments; i++) {
            const angle = i / segments * 2 * Math.PI;
            // Full revolution
            const point = this.profileCurve.evaluate(angle);
            // Evaluate the profile curve
            controlPoints.push([
                point[0] * Math.cos(angle),
                // X position after rotation
                point[0] * Math.sin(angle),
                // Y position after rotation
                point[1]
            ]);
        }
        // Z position remains the same
        this.controlPoints = controlPoints;
    }
    evaluate(u, v) {
        const point = this.profileCurve.evaluate(u * 2 * Math.PI);
        // Evaluate with angle
        const angle = v * 2 * Math.PI;
        // Angle for the revolution
        return [
            point[0] * Math.cos(angle),
            point[0] * Math.sin(angle),
            point[1]
        ];
    }
}
class SubdivisionSurface extends NURBS_Surface {
    constructor(controlPoints) {
        super(2, 2, controlPoints, [], []);
    }
    evaluate(u, v) {
        const controlPoints = this.controlPoints;
        const point = [
            0,
            0,
            0
        ];
        const degreeU = this.degreeU;
        const degreeV = this.degreeV;
        for (let i = 0; i <= degreeU; i++) {
            for (let j = 0; j <= degreeV; j++) {
                const basisU = this.basisFunctionU(i, degreeU, u);
                const basisV = this.basisFunctionV(j, degreeV, v);
                const controlPoint = controlPoints[i][j];
                for (let k = 0; k < point.length; k++) {
                    point[k] += basisU * basisV * controlPoint[k];
                }
            }
        }
        return point;
    }
    getPoints(stepSizeU, stepSizeV) {
        const points = [];
        for (let u = 0; u <= 1; u += stepSizeU) {
            for (let v = 0; v <= 1; v += stepSizeV) {
                points.push(this.evaluate(u, v));
            }
        }
        return points;
    }
    generateControlPoints() {
        const controlPoints = [];
        for (let i = 0; i <= 1; i++) {
            const row = [];
            for (let j = 0; j <= 1; j++) {
                row.push([
                    i,
                    j,
                    0
                ]);
            }
            // Example logic for generating control points
            controlPoints.push(row);
        }
        this.controlPoints = controlPoints;
    }
}
class VectorUtils {
    static addVectors(v1, v2) {
        return v1.map((val, index) => val + v2[index]);
    }
    static subtractVectors(v1, v2) {
        return v1.map((val, index) => val - v2[index]);
    }
    static normalizeVector(v) {
        const length = Math.sqrt(v.reduce((sum, val) => sum + val ** 2, 0));
        if (length === 0) {
            throw new Error('Cannot normalize a zero-length vector.');
        }
        return v.map(val => val / length);
    }
    static isClose(point1, point2, tolerance = 0.000001) {
        return point1.every((val, index) => Math.abs(val - point2[index]) < tolerance);
    }
}
class BREP {
    constructor() {
        this.faces = [];
    }
    addFace(nurbsSurface) {
        const face = new Face(nurbsSurface);
        this.faces.push(face);
        return face;
    }
    getFaces() {
        return this.faces;
    }
    getEdges() {
        return this.faces.flatMap(face => face.getEdges());
    }
    getVertices() {
        const verticesSet = new Set();
        this.faces.forEach(face => {
            face.getVertices().forEach(vertex => {
                verticesSet.add(vertex.getCoordinates().toString());
            });
        });
        return Array.from(verticesSet).map(vertexString => vertexString.split(',').map(Number));
    }
    splitEdge(edge, parameter) {
        const controlPoints = edge.nurbsCurve.controlPoints;
        const splitPoint = edge.evaluate(parameter);
        const leftControlPoints = controlPoints.slice(0, Math.floor(controlPoints.length / 2)).concat(splitPoint);
        const rightControlPoints = [splitPoint].concat(controlPoints.slice(Math.floor(controlPoints.length / 2)));
        const degree = edge.nurbsCurve.degree;
        const leftCurve = new NURBS_Curve(degree, leftControlPoints, edge.nurbsCurve.knots);
        const rightCurve = new NURBS_Curve(degree, rightControlPoints, edge.nurbsCurve.knots);
        return {
            leftEdge: new Edge(leftCurve),
            rightEdge: new Edge(rightCurve)
        };
    }
    removeFace(face) {
        const index = this.faces.indexOf(face);
        if (index !== -1) {
            this.faces.splice(index, 1);
        } else {
            throw new Error('Face not found in the BREP.');
        }
    }
    mergeFaces(face1, face2) {
        // Placeholder logic for merging faces
        return new Face(face1.nurbsSurface);
    }
    splitFaceByPlane(face, point, normal) {
        const edges = face.getEdges();
        const aboveFace = new Face(face.nurbsSurface);
        const belowFace = new Face(face.nurbsSurface);
        edges.forEach(edge => {
            const midPoint = edge.getMidPoint();
            if (this.isPointAbovePlane(midPoint, point, normal)) {
                aboveFace.addEdge(edge);
            } else {
                belowFace.addEdge(edge);
            }
        });
        this.removeFace(face);
        this.addFace(aboveFace);
        this.addFace(belowFace);
    }
    isPointAbovePlane(point, planePoint, normal) {
        const vectorToPoint = VectorUtils.subtractVectors(point, planePoint);
        const dotProduct = vectorToPoint.reduce((sum, val, index) => sum + val * normal[index], 0);
        return dotProduct > 0;
    }
    getFaceNormals() {
        return this.faces.map(face => face.getNormal());
    }
    removeVertexFromEdges(vertex) {
        this.faces.forEach(face => {
            face.getEdges().forEach(edge => {
                edge.nurbsCurve.controlPoints = edge.nurbsCurve.controlPoints.filter(cp => !this.isClose(cp, vertex.getCoordinates()));
            });
        });
    }
    mergeEdgesByVertices(vertex1, vertex2) {
        const edge = this.findEdgeWithVertices(vertex1, vertex2);
        if (!edge) {
            throw new Error('Cannot merge edges that do not share vertices.');
        }
        const combinedControlPoints = [
            ...edge.nurbsCurve.controlPoints.slice(0, -1),
            ...edge.nurbsCurve.controlPoints
        ];
        const degree = Math.max(edge.nurbsCurve.degree, edge.nurbsCurve.degree);
        const newKnots = this.mergeKnotVectors(edge.nurbsCurve.knots, edge.nurbsCurve.knots);
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    combineFaces(face1, face2) {
        if (!this.isFaceShared(face1, face2)) {
            throw new Error('Faces do not share edges and cannot be combined.');
        }
        const combinedControlPoints = [
            ...face1.nurbsSurface.controlPoints,
            ...face2.nurbsSurface.controlPoints
        ];
        const degreeU = Math.max(face1.nurbsSurface.degreeU, face2.nurbsSurface.degreeU);
        const degreeV = Math.max(face1.nurbsSurface.degreeV, face2.nurbsSurface.degreeV);
        const newKnotsU = this.mergeKnotVectors(face1.nurbsSurface.knotsU, face2.nurbsSurface.knotsU);
        const newKnotsV = this.mergeKnotVectors(face1.nurbsSurface.knotsV, face2.nurbsSurface.knotsV);
        const combinedSurface = new NURBS_Surface(degreeU, degreeV, combinedControlPoints, newKnotsU, newKnotsV);
        const combinedFace = new Face(combinedSurface);
        this.removeFace(face1);
        this.removeFace(face2);
        this.addFace(combinedFace);
        return combinedFace;
    }
    addVertexToFace(face, vertex) {
        const edge = new Edge(vertex.getCoordinates());
        face.addEdge(edge);
        return vertex;
    }
    splitEdgeByPoint(edge, point) {
        const curve = edge.nurbsCurve;
        const t = this.findParameterAtPoint(curve, point);
        if (t === null) {
            throw new Error('Point is not on the curve.');
        }
        const splitCurve = this.splitEdge(edge, t);
        edge.nurbsCurve = splitCurve.leftCurve;
        // Keep the left part, for example
        this.addEdge(splitCurve.rightCurve);
    }
    findParameterAtPoint(curve, point) {
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const curvePoint = curve.evaluate(t);
            if (this.isClose(curvePoint, point)) {
                return t;
            }
        }
        return null;
    }
    isClose(point1, point2, tolerance = 0.000001) {
        return point1.every((val, index) => Math.abs(val - point2[index]) < tolerance);
    }
    findEdgeWithVertices(vertex1, vertex2) {
        const edges = this.getEdges();
        return edges.find(edge => {
            const controlPoints = edge.nurbsCurve.controlPoints;
            return this.isClose(controlPoints[0], vertex1.getCoordinates()) && this.isClose(controlPoints[1], vertex2.getCoordinates()) || this.isClose(controlPoints[0], vertex2.getCoordinates()) && this.isClose(controlPoints[1], vertex1.getCoordinates());
        }) || null;
    }
    splitEdgeByParameter(edge, parameter) {
        const splitResult = edge.splitAtParameter(parameter);
        this.removeEdge(edge);
        return [
            splitResult.leftEdge,
            splitResult.rightEdge
        ];
    }
    getFaceByIndex(index) {
        if (index < 0 || index >= this.faces.length) {
            throw new Error('Face index is out of bounds.');
        }
        return this.faces[index];
    }
    findFaceBySurface(nurbsSurface) {
        return this.faces.find(face => face.nurbsSurface === nurbsSurface) || null;
    }
    mergeAdjacentFaces() {
        const mergedFaces = [];
        for (let i = 0; i < this.faces.length; i++) {
            for (let j = i + 1; j < this.faces.length; j++) {
                const sharedEdges = this.findSharedEdges(this.faces[i], this.faces[j]);
                if (sharedEdges.length > 0) {
                    const mergedFace = this.mergeFaces(this.faces[i], this.faces[j]);
                    mergedFaces.push(mergedFace);
                    this.removeFace(this.faces[j]);
                    this.removeFace(this.faces[i]);
                    i--;
                    break;
                }
            }
        }
        mergedFaces.forEach(face => this.addFace(face));
    }
    splitFace(face, curve) {
        const edges = face.getEdges();
        const splitEdges = edges.filter(edge => curve.isClose(edge.nurbsCurve));
        const newControlPoints1 = [];
        const newControlPoints2 = [];
        splitEdges.forEach(edge => {
            const controlPoints = edge.nurbsCurve.controlPoints;
            const midPointIndex = Math.floor(controlPoints.length / 2);
            newControlPoints1.push(...controlPoints.slice(0, midPointIndex));
            newControlPoints2.push(...controlPoints.slice(midPointIndex));
        });
        const newSurface1 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints1, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newSurface2 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints2, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newFace1 = new Face(newSurface1);
        const newFace2 = new Face(newSurface2);
        this.removeFace(face);
        this.addFace(newFace1);
        this.addFace(newFace2);
        return [
            newFace1,
            newFace2
        ];
    }
    getEdgesByFace(face) {
        return face.getEdges();
    }
    mergeEdgesByFaces(face1, face2) {
        const edgesFace1 = this.getEdgesByFace(face1);
        const edgesFace2 = this.getEdgesByFace(face2);
        const sharedEdges = edgesFace1.filter(edge1 => edgesFace2.some(edge2 => this.isEdgeShared(edge1, edge2)));
        if (sharedEdges.length === 0) {
            throw new Error('No shared edges between the faces, cannot merge.');
        }
        const combinedControlPoints = sharedEdges.reduce((points, edge) => {
            return [
                ...points,
                ...edge.nurbsCurve.controlPoints
            ];
        }, []);
        const degree = Math.max(...sharedEdges.map(edge => edge.nurbsCurve.degree));
        const newKnots = this.mergeKnotVectors(sharedEdges.map(edge => edge.nurbsCurve.knots));
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    findFacesByEdge(edge) {
        return this.faces.filter(face => face.getEdges().includes(edge));
    }
    isSurfaceClosed(surface) {
        const borderEdges = surface.getEdges();
        return borderEdges.every(edge => edge.getControlPoints()[0] === edge.getControlPoints()[1]);
    }
    updateVertex(vertex, newCoordinates) {
        const foundVertex = this.findVertex(vertex);
        if (foundVertex) {
            foundVertex.point = newCoordinates;
            this.getEdgesByVertex(foundVertex).forEach(edge => {
                edge.nurbsCurve.controlPoints = edge.nurbsCurve.controlPoints.map(cp => this.isClose(cp, vertex.getCoordinates()) ? newCoordinates : cp);
            });
        } else {
            throw new Error('Vertex not found.');
        }
    }
    updateConnectedEdges(vertex) {
        const edges = this.getEdgesByVertex(vertex);
        edges.forEach(edge => {
            const curve = edge.nurbsCurve;
            curve.controlPoints = curve.controlPoints.map(cp => {
                if (this.isClose(cp, vertex.getCoordinates())) {
                    return vertex.getCoordinates();
                }
                return cp;
            });
        });
    }
    findVertex(vertex) {
        for (const face of this.faces) {
            const vertices = face.getVertices();
            const foundVertex = vertices.find(v => this.isClose(v.getCoordinates(), vertex.getCoordinates()));
            if (foundVertex) {
                return foundVertex;
            }
        }
        return null;
    }
    getEdgesByVertex(vertex) {
        const edges = [];
        for (const face of this.faces) {
            for (const edge of face.getEdges()) {
                if (edge.nurbsCurve.controlPoints.some(cp => this.isClose(cp, vertex.getCoordinates()))) {
                    edges.push(edge);
                }
            }
        }
        return edges;
    }
    clearFaces() {
        this.faces = [];
    }
    updateFace(face, newNurbsSurface) {
        face.nurbsSurface = newNurbsSurface;
    }
    getAllVertices() {
        const vertices = new Set();
        this.faces.forEach(face => {
            face.getVertices().forEach(vertex => {
                vertices.add(vertex.getCoordinates().toString());
            });
        });
        return Array.from(vertices).map(vertexString => vertexString.split(',').map(Number));
    }
    isValid() {
        return this.faces.length > 0 && this.faces.every(face => face.getEdges().length > 0);
    }
    combineEdges(edgeList) {
        if (!edgeList || edgeList.length < 2) {
            throw new Error('At least two edges are required for combination.');
        }
        const combinedControlPoints = edgeList.flatMap(edge => edge.nurbsCurve.controlPoints);
        const degree = Math.max(...edgeList.map(edge => edge.nurbsCurve.degree));
        const newKnots = this.mergeKnotVectors(edgeList.map(edge => edge.nurbsCurve.knots));
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    edgesShareVertices(edge1, edge2) {
        const verticesEdge1 = edge1.nurbsCurve.controlPoints;
        const verticesEdge2 = edge2.nurbsCurve.controlPoints;
        return verticesEdge1.some(v1 => verticesEdge2.some(v2 => this.isClose(v1, v2)));
    }
    getFaceByVertices(vertex1, vertex2) {
        for (const face of this.faces) {
            const vertices = face.getVertices();
            if (vertices.some(v => this.isClose(v.getCoordinates(), vertex1)) && vertices.some(v => this.isClose(v.getCoordinates(), vertex2))) {
                return face;
            }
        }
        return null;
    }
    mergeFacesByVertices(vertex1, vertex2) {
        const face1 = this.getFaceByVertices(vertex1, vertex2);
        const face2 = this.getFaceByVertices(vertex2, vertex1);
        if (!face1 || !face2) {
            throw new Error('No faces are sharing the specified vertices.');
        }
        return this.mergeFaces(face1, face2);
    }
    findVertexByCoordinates(coordinates) {
        for (const face of this.faces) {
            const vertices = face.getVertices();
            const foundVertex = vertices.find(vertex => this.isClose(vertex.getCoordinates(), coordinates));
            if (foundVertex) {
                return foundVertex;
            }
        }
        return null;
    }
    getEdgesSharingVertex(vertex) {
        const edges = [];
        for (const face of this.faces) {
            const faceEdges = face.getEdges();
            for (const edge of faceEdges) {
                if (edge.getControlPoints().some(cp => this.isClose(cp, vertex.getCoordinates()))) {
                    edges.push(edge);
                }
            }
        }
        return edges;
    }
    splitEdgeByVertices(vertex1, vertex2) {
        const edge = this.findEdgeWithVertices(vertex1, vertex2);
        if (!edge) {
            throw new Error('No edge found between the specified vertices.');
        }
        const midPoint = edge.getMidPoint();
        const newCurve = this.splitEdge(edge, 0.5);
        return newCurve;
    }
    getFaceEdges(face) {
        return face.getEdges();
    }
    mergeSharedFaces() {
        const mergedFaces = [];
        for (let i = 0; i < this.faces.length; i++) {
            for (let j = i + 1; j < this.faces.length; j++) {
                const sharedEdges = this.findSharedEdges(this.faces[i], this.faces[j]);
                if (sharedEdges.length > 0) {
                    const mergedFace = this.mergeFaces(this.faces[i], this.faces[j]);
                    mergedFaces.push(mergedFace);
                    this.removeFace(this.faces[i]);
                    this.removeFace(this.faces[j]);
                    i--;
                    break;
                }
            }
        }
        mergedFaces.forEach(face => this.addFace(face));
    }
    splitFaceByCurve(face, curve) {
        const newEdges = [];
        const intersections = this.intersectEdges(face.getEdges(), curve);
        intersections.forEach(intersection => {
            const splitResult = this.splitEdgeByParameter(intersection.edge, intersection.param1);
            newEdges.push(splitResult.leftEdge);
            newEdges.push(splitResult.rightEdge);
        });
        const newFace = new Face(face.nurbsSurface);
        newEdges.forEach(edge => newFace.addEdge(edge));
        this.removeFace(face);
        this.addFace(newFace);
        return newFace;
    }
    findSharedEdges(face1, face2) {
        return face1.getEdges().filter(edge1 => face2.getEdges().some(edge2 => this.isEdgeShared(edge1, edge2)));
    }
    splitEdgeByRatio(edge, ratio) {
        const controlPoints = edge.nurbsCurve.controlPoints;
        const midIndex = Math.floor(controlPoints.length * ratio);
        const firstHalf = controlPoints.slice(0, midIndex);
        const secondHalf = controlPoints.slice(midIndex);
        edge.setControlPoints(firstHalf);
        const newEdge = new Edge(new NURBS_Curve(edge.nurbsCurve.degree, secondHalf, edge.nurbsCurve.knots));
        this.addEdge(newEdge);
    }
    addEdgeToFace(face, nurbsCurve) {
        const edge = new Edge(nurbsCurve);
        face.addEdge(edge);
        return edge;
    }
    weldEdges(edge1, edge2) {
        if (!this.edgesShareVertices(edge1, edge2)) {
            throw new Error('Cannot weld edges that do not share vertices.');
        }
        const combinedControlPoints = [
            ...edge1.nurbsCurve.controlPoints.slice(0, -1),
            ...edge2.nurbsCurve.controlPoints
        ];
        const degree = Math.max(edge1.nurbsCurve.degree, edge2.nurbsCurve.degree);
        const newKnots = this.mergeKnotVectors(edge1.nurbsCurve.knots, edge2.nurbsCurve.knots);
        const weldedCurve = new NURBS_Curve(degree, combinedControlPoints, newKnots);
        return new Edge(weldedCurve);
    }
    getTopologicalInfo() {
        return {
            vertexCount: this.getVertexCount(),
            edgeCount: this.getEdgeCount(),
            faceCount: this.getFaceCount()
        };
    }
    containsPoint(point) {
        return this.BREP.getFaces().some(face => {
            const vertices = face.getVertices();
            return vertices.some(vertex => this.isClose(vertex.getCoordinates(), point));
        });
    }
    removeEdgeFromFace(face, edge) {
        const edges = face.getEdges();
        const index = edges.indexOf(edge);
        if (index !== -1) {
            edges.splice(index, 1);
        } else {
            throw new Error('Edge does not exist in the face.');
        }
    }
    isEdgeShared(edge1, edge2) {
        const verticesEdge1 = edge1.getVertices().map(v => v.getCoordinates());
        const verticesEdge2 = edge2.getVertices().map(v => v.getCoordinates());
        return verticesEdge1.some(v1 => verticesEdge2.some(v2 => this.isClose(v1, v2)));
    }
    splitFaceByEdge(face, edge) {
        const edges = face.getEdges();
        const splitEdges = edges.filter(e => e === edge);
        const newControlPoints1 = [];
        const newControlPoints2 = [];
        splitEdges.forEach(edge => {
            const controlPoints = edge.nurbsCurve.controlPoints;
            const midPointIndex = Math.floor(controlPoints.length / 2);
            newControlPoints1.push(...controlPoints.slice(0, midPointIndex + 1));
            newControlPoints2.push(...controlPoints.slice(midPointIndex));
        });
        const newSurface1 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints1, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newSurface2 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints2, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newFace1 = new Face(newSurface1);
        const newFace2 = new Face(newSurface2);
        this.removeFace(face);
        this.addFace(newFace1);
        this.addFace(newFace2);
        return [
            newFace1,
            newFace2
        ];
    }
    getEdgeCount() {
        return this.faces.reduce((count, face) => count + face.getEdges().length, 0);
    }
    getVertexCount() {
        return this.faces.reduce((count, face) => count + face.getVertices().length, 0);
    }
    getFaceCount() {
        return this.faces.length;
    }
    mergeKnotVectors(knots1, knots2) {
        const combinedKnots = new Set([
            ...knots1,
            ...knots2
        ]);
        return [...combinedKnots].sort((a, b) => a - b);
    }
    clear() {
        this.faces.forEach(face => face.clearEdges());
        this.faces = [];
    }
    isFaceValid(face) {
        return face.getEdges().length > 0 && face.getEdges().every(edge => edge.isValid());
    }
    areFacesConnected(face1, face2) {
        return this.findSharedEdges(face1, face2).length > 0;
    }
    isFaceShared(face1, face2) {
        const sharedEdges = this.findSharedEdges(face1, face2);
        return sharedEdges.length > 0;
    }
    getEdgeVertices(edge) {
        return edge.nurbsCurve.controlPoints.map(cp => new Vertex(cp));
    }
    mergeEdges(edge1, edge2) {
        if (!this.edgesShareVertices(edge1, edge2)) {
            throw new Error('Cannot merge edges that do not share vertices.');
        }
        const combinedControlPoints = [
            ...edge1.nurbsCurve.controlPoints.slice(0, -1),
            ...edge2.nurbsCurve.controlPoints
        ];
        const degree = Math.max(edge1.nurbsCurve.degree, edge2.nurbsCurve.degree);
        const newKnots = this.mergeKnotVectors(edge1.nurbsCurve.knots, edge2.nurbsCurve.knots);
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    weldEdge(edge) {
        const adjacentFaces = this.findFacesByEdge(edge);
        if (adjacentFaces.length !== 2) {
            throw new Error('Edge must be shared by exactly two faces to be welded.');
        }
        const mergedFace = this.mergeFaces(adjacentFaces[0], adjacentFaces[1]);
        this.removeFace(adjacentFaces[0]);
        this.removeFace(adjacentFaces[1]);
        this.addFace(mergedFace);
        return mergedFace;
    }
    removeOrphanedEdges() {
        const allEdges = this.getEdges();
        const usedEdges = new Set(this.getFaces().flatMap(f => f.getEdges()));
        const orphanedEdges = allEdges.filter(edge => !usedEdges.has(edge));
        orphanedEdges.forEach(orphanedEdge => {
            this.faces.forEach(face => {
                face.removeEdge(orphanedEdge);
            });
        });
    }
    isFaceConnected(face1, face2) {
        return this.findSharedEdges(face1, face2).length > 0;
    }
    findOrCreateVertex(point) {
        const existingVertex = this.findVertexByCoordinates(point);
        if (existingVertex) {
            return existingVertex;
        }
        const newVertex = new Vertex(point);
        this.addVertexToFace(newFace, newVertex);
        return newVertex;
    }
    calculateVolume() {
        let totalVolume = 0;
        for (const face of this.faces) {
            const area = face.getArea();
            const normal = face.getNormal();
            const referencePoint = face.nurbsSurface.controlPoints[0][0];
            // arbitrary reference point
            const dotProduct = normal.reduce((sum, val, index) => sum + val * referencePoint[index], 0);
            const signedVolume = 1 / 3 * area * dotProduct;
            totalVolume += Math.abs(signedVolume);
        }
        return totalVolume;
    }
    computeVolume() {
        return this.faces.reduce((totalVolume, face) => totalVolume + face.getArea() * face.height, 0);
    }
    computeSurfaceArea() {
        return this.faces.reduce((totalArea, face) => totalArea + face.getArea(), 0);
    }
    isFaceVisible(face) {
        const normal = face.getNormal();
        const viewingDirection = [
            0,
            0,
            1
        ];
        // Assuming viewing direction along the positive z-axis
        const dotProduct = normal.reduce((sum, val, index) => sum + val * viewingDirection[index], 0);
        return dotProduct > 0;
    }
    calculateSurfaceArea() {
        return this.faces.reduce((totalArea, face) => totalArea + face.getArea(), 0);
    }
    getEdgeByVertices(vertex1, vertex2) {
        for (const face of this.faces) {
            const edges = face.getEdges();
            for (const edge of edges) {
                if (this.isClose(edge.nurbsCurve.controlPoints[0], vertex1.getCoordinates()) && this.isClose(edge.nurbsCurve.controlPoints[1], vertex2.getCoordinates()) || this.isClose(edge.nurbsCurve.controlPoints[0], vertex2.getCoordinates()) && this.isClose(edge.nurbsCurve.controlPoints[1], vertex1.getCoordinates())) {
                    return edge;
                }
            }
        }
        return null;
    }
    getAdjacentFaces(face) {
        return this.faces.filter(f => f.getEdges().some(edge => face.getEdges().includes(edge)));
    }
    validate() {
        const faceValidity = this.faces.every(face => face.isFaceValid());
        const edgeValidity = this.getEdges().every(edge => edge.isValid());
        return faceValidity && edgeValidity;
    }
    isClosed() {
        return this.faces.length > 0 && this.faces.every(face => face.getEdges().every(edge => edge.nurbsCurve.controlPoints[0] === edge.nurbsCurve.controlPoints[1]));
    }
    updateEdge(edge, newNurbsCurve) {
        const existingEdge = this.getEdges().find(e => e === edge);
        if (existingEdge) {
            existingEdge.nurbsCurve = newNurbsCurve;
        } else {
            throw new Error('Edge not found for update.');
        }
    }
    clearBREP() {
        this.faces.forEach(face => face.clearEdges());
        this.faces = [];
    }
    computeFaceNormals() {
        this.faces.forEach(face => {
            const normal = face.calculateNormal();
            face.normal = normal;
        });
    }
    isEdgeValid(edge) {
        // Validates an edge by checking if it exists and has more than one control point
        return edge && edge.nurbsCurve.controlPoints.length > 1;
    }
    mergeFacesByEdge(face1, face2) {
        const sharedEdges = this.findSharedEdges(face1, face2);
        if (sharedEdges.length === 0) {
            throw new Error('No shared edges between the faces, cannot merge.');
        }
        const combinedControlPoints = sharedEdges.reduce((points, edge) => {
            return [
                ...points,
                ...edge.nurbsCurve.controlPoints
            ];
        }, []);
        const degree = Math.max(face1.nurbsSurface.degreeU, face2.nurbsSurface.degreeU);
        const newKnotsU = this.mergeKnotVectors(face1.nurbsSurface.knotsU, face2.nurbsSurface.knotsU);
        const newKnotsV = this.mergeKnotVectors(face1.nurbsSurface.knotsV, face2.nurbsSurface.knotsV);
        const mergedSurface = new NURBS_Surface(degree, degree, combinedControlPoints, newKnotsU, newKnotsV);
        return new Face(mergedSurface);
    }
    combineCurves(curve1, curve2) {
        const combinedControlPoints = [
            ...curve1.controlPoints,
            ...curve2.controlPoints.slice(1)
        ];
        // Avoid duplicate endpoint
        const degree = Math.max(curve1.degree, curve2.degree);
        const newKnots = this.mergeKnotVectors(curve1.knots, curve2.knots);
        return new NURBS_Curve(degree, combinedControlPoints, newKnots);
    }
    getFacesByEdge(edge) {
        return this.faces.filter(face => face.getEdges().includes(edge));
    }
    removeVertex(vertex) {
        this.faces.forEach(face => {
            face.getEdges().forEach(edge => {
                edge.nurbsCurve.controlPoints = edge.nurbsCurve.controlPoints.filter(cp => !this.isClose(cp, vertex.getCoordinates()));
            });
        });
    }
    removeEdge(edge) {
        for (const face of this.faces) {
            const edges = face.getEdges();
            const index = edges.indexOf(edge);
            if (index !== -1) {
                edges.splice(index, 1);
                return;
            }
        }
        throw new Error('Edge not found for removal.');
    }
    findFaceByEdge(edge) {
        // Finds the face that contains the given edge
        return this.faces.find(face => face.getEdges().includes(edge)) || null;
    }
    // Replace with actual merging logic
    addEdge(edge) {
    }
    union(brep1, brep2) {
        if (!brep1 || !brep2) {
            throw new Error('Both BREP instances are required');
        }
        const unionBREP = new BREP();
        // Add faces from both BREP instances
        brep1.getFaces().forEach(face => unionBREP.addFace(face));
        brep2.getFaces().forEach(face => unionBREP.addFace(face));
        // Remove intersecting faces
        const intersectingFaces = this.identifyIntersections(brep1, brep2);
        intersectingFaces.forEach(({face1}) => {
            unionBREP.removeFace(face1);
        });
        return unionBREP;
    }
    intersection(brep1, brep2) {
        const intersections = this.identifyIntersections(brep1, brep2);
        const intersectionBREP = new BREP();
        intersections.forEach(({face1, face2}) => {
            const sharedEdges = this.findSharedEdges(face1, face2);
            if (sharedEdges.length > 0) {
                const intersectionSurface = this.createIntersectionSurface(face1, face2, sharedEdges);
                if (intersectionSurface) {
                    intersectionBREP.addFace(intersectionSurface);
                    sharedEdges.forEach(edge => {
                        const newEdgeCurve = this.processIntersectedEdge(edge, intersectionSurface);
                        const newEdge = new Edge(newEdgeCurve);
                        if (intersectionBREP.getFaces().length > 0) {
                            intersectionBREP.getFaces().at(-1).addEdge(newEdge);
                        }
                    });
                }
            }
        });
        return intersectionBREP;
    }
    difference(brep1, brep2) {
        if (!brep1 || !brep2) {
            throw new Error('Both BREP instances are required');
        }
        const intersectionBREP = this.intersection(brep1, brep2);
        const differenceBREP = new BREP();
        brep1.getFaces().forEach(face => {
            if (!intersectionBREP.containsFace(face)) {
                differenceBREP.addFace(face);
            }
        });
        return differenceBREP;
    }
    identifyIntersections(brep1, brep2) {
        const intersections = [];
        for (const face1 of brep1.getFaces()) {
            for (const face2 of brep2.getFaces()) {
                const faceIntersections = this.intersectFaces(face1, face2);
                if (faceIntersections.length > 0) {
                    intersections.push({
                        face1,
                        face2,
                        intersections: faceIntersections
                    });
                }
            }
        }
        return intersections;
    }
    modifyTopology() {
        // Initialize variables to keep track of changes
        const modifiedEdges = [];
        const modifiedFaces = [];
        // Iterate through the faces
        this.faces.forEach(face => {
            const edges = face.getEdges();
            const sharedEdges = edges.filter(edge => this.findSharedFacesByEdge(edge).length > 1);
            // If shared edges are found, these might need to be modified
            if (sharedEdges.length > 0) {
                sharedEdges.forEach(edge => {
                    // Split the shared edge
                    this.splitEdgeByRatio(edge, 0.5);
                    // Keep track of the modified edges
                    modifiedEdges.push(edge);
                    // Find the adjacent faces that share this edge
                    const adjacentFaces = this.findFacesByEdge(edge);
                    if (adjacentFaces.length > 1) {
                        // Modify the topology by merging faces if they can be combined
                        const mergedFace = this.mergeFaces(adjacentFaces[0], adjacentFaces[1]);
                        modifiedFaces.push(mergedFace);
                    }
                });
            }
        });
        // Return the modified edges and faces for further processing or review
        return {
            modifiedEdges,
            modifiedFaces
        };
    }
    createIntersectionSurface(face1, face2, sharedEdges) {
        const combinedControlPoints = sharedEdges.reduce((acc, edge) => acc.concat(edge.nurbsCurve.controlPoints), []);
        const degreeU = Math.max(face1.nurbsSurface.degreeU, face2.nurbsSurface.degreeU);
        const degreeV = Math.max(face1.nurbsSurface.degreeV, face2.nurbsSurface.degreeV);
        const newKnotsU = this.mergeKnotVectors(face1.nurbsSurface.knotsU, face2.nurbsSurface.knotsU);
        const newKnotsV = this.mergeKnotVectors(face1.nurbsSurface.knotsV, face2.nurbsSurface.knotsV);
        return new Face(new NURBS_Surface(degreeU, degreeV, combinedControlPoints, newKnotsU, newKnotsV));
    }
    processIntersectedEdge(edge, intersectionSurface) {
        const span = edge.nurbsCurve.findSpan(edge.nurbsCurve.degree, 0.5, edge.nurbsCurve.knots);
        const newControlPoints = edge.nurbsCurve.controlPoints.map((cp, index) => {
            const basePoint = intersectionSurface.nurbsSurface.evaluate(cp[0], cp[1]);
            return edge.nurbsCurve.blendControlPoints(cp, basePoint, 0.5, span, edge.nurbsCurve.degree);
        });
        return new NURBS_Curve(edge.nurbsCurve.degree, newControlPoints, edge.nurbsCurve.knots);
    }
    splitEdgeByEdge(edge, intersectingEdge) {
        // Split the edge at all intersection points with the intersectingEdge
        const splitEdges = [];
        const intersectionPoints = this.findEdgeIntersections(edge, intersectingEdge);
        let currentEdge = edge;
        intersectionPoints.forEach(point => {
            const parameter = this.findParameterAtPoint(currentEdge.nurbsCurve, point);
            if (parameter !== null) {
                const splitResult = currentEdge.splitAtParameter(parameter);
                splitEdges.push(splitResult.leftEdge);
                currentEdge = splitResult.rightEdge;
            }
        });
        splitEdges.push(currentEdge);
        return splitEdges;
    }
    findEdgeIntersections(edge1, edge2) {
        const intersections = this.intersectEdges(edge1, edge2);
        if (intersections) {
            intersections.forEach(intersection => {
                const edge1Param = this.findEdgeParameter(edge1, intersection.point);
                const edge2Param = this.findEdgeParameter(edge2, intersection.point);
                intersection.edge1Param = edge1Param;
                intersection.edge2Param = edge2Param;
            });
        }
        return intersections || [];
    }
    removeIntersectingFaces(baseFace) {
        const intersectingFaces = this.faces.filter(face => {
            if (face === baseFace)
                return false;
            return this.intersectFaces(baseFace, face).length > 0;
        });
        intersectingFaces.forEach(face => this.removeFace(face));
    }
    containsFace(face) {
        return this.faces.some(f => f.nurbsSurface === face.nurbsSurface);
    }
    isFaceIdentical(face1, face2) {
        const face1Edges = face1.getEdges().map(edge => edge.nurbsCurve);
        const face2Edges = face2.getEdges().map(edge => edge.nurbsCurve);
        if (face1Edges.length !== face2Edges.length) {
            return false;
        }
        return face1Edges.every(curve1 => face2Edges.some(curve2 => this.isCurveIdentical(curve1, curve2)));
    }
    isCurveIdentical(curve1, curve2) {
        return curve1.degree === curve2.degree && this.areControlPointsIdentical(curve1.controlPoints, curve2.controlPoints) && this.areKnotsIdentical(curve1.knots, curve2.knots);
    }
    areControlPointsIdentical(cp1, cp2) {
        if (cp1.length !== cp2.length) {
            return false;
        }
        return cp1.every((point1, index) => this.isClose(point1, cp2[index]));
    }
    areKnotsIdentical(knots1, knots2) {
        if (knots1.length !== knots2.length) {
            return false;
        }
        return knots1.every((knot, index) => knot === knots2[index]);
    }
    findSharedFacesByEdge(edge) {
        return this.faces.filter(face => face.getEdges().includes(edge));
    }
    checkIntersection(face1, face2) {
        return face1.getEdges().some(edge1 => face2.getEdges().some(edge2 => this.intersectEdges(edge1, edge2).length > 0));
    }
    subtractBREPs(brep1, brep2) {
        // Check for valid input
        if (!brep1 || !brep2) {
            throw new Error('Both BREP instances are required');
        }
        // Find intersecting faces between brep1 and brep2
        const intersectionBREP = this.intersection(brep1, brep2);
        // Initialize the result BREP
        const resultBREP = new BREP();
        // Add all faces from brep1 that do not intersect with brep2
        brep1.getFaces().forEach(face => {
            if (!intersectionBREP.containsFace(face)) {
                resultBREP.addFace(face);
            }
        });
        // Logically handle potential remaining surface differences here
        // Using remaining intersections or offset surfaces can be implemented.
        return resultBREP;
    }
    subtractBREPS(brep1, brep2) {
        if (!brep1 || !brep2) {
            throw new Error('Both BREP instances are required');
        }
        const intersectionBREP = this.intersection(brep1, brep2);
        const resultBREP = new BREP();
        brep1.getFaces().forEach(face => {
            if (!intersectionBREP.containsFace(face)) {
                resultBREP.addFace(face);
            }
        });
        return resultBREP;
    }
    intersectFaceWithPlane(face, planePoint, planeNormal) {
        const intersections = [];
        face.getEdges().forEach(edge => {
            const intersectionPoints = this.intersectEdgeWithPlane(edge, planePoint, planeNormal);
            if (intersectionPoints) {
                intersections.push(...intersectionPoints.map(point => ({ point })));
            }
        });
        return intersections;
    }
    intersectEdgeWithPlane(edge, planePoint, planeNormal) {
        const intersections = [];
        const controlPoints = edge.getControlPoints();
        for (let i = 0; i < controlPoints.length - 1; i++) {
            const p0 = controlPoints[i];
            const p1 = controlPoints[i + 1];
            const direction = VectorUtils.subtractVectors(p1, p0);
            const dotProductDirectionNormal = planeNormal.reduce((sum, val, index) => sum + val * direction[index], 0);
            if (Math.abs(dotProductDirectionNormal) > Number.EPSILON) {
                const t = planeNormal.reduce((sum, val, index) => sum + val * (planePoint[index] - p0[index]), 0) / dotProductDirectionNormal;
                if (t >= 0 && t <= 1) {
                    intersections.push(p0.map((coord, index) => coord + t * direction[index]));
                }
            }
        }
        return intersections.length > 0 ? intersections : null;
    }
    areEdgesIntersecting(edge1, edge2) {
        return this.intersectEdges(edge1, edge2).length > 0;
    }
    intersectEdges(edge1, edge2) {
        const intersections = [];
        const stepSize = 0.01;
        for (let t1 = 0; t1 <= 1; t1 += stepSize) {
            const point1 = edge1.evaluate(t1);
            for (let t2 = 0; t2 <= 1; t2 += stepSize) {
                const point2 = edge2.evaluate(t2);
                if (VectorUtils.isClose(point1, point2)) {
                    intersections.push({
                        point: point1,
                        param1: t1,
                        param2: t2
                    });
                }
            }
        }
        return intersections;
    }
    intersectFaces(face1, face2) {
        const intersections = [];
        face1.getEdges().forEach(edge1 => {
            face2.getEdges().forEach(edge2 => {
                const edgeIntersections = edge1.intersectEdges(edge2);
                if (edgeIntersections.length > 0) {
                    intersections.push({
                        edge1,
                        edge2,
                        points: edgeIntersections
                    });
                }
            });
        });
        return intersections;
    }
    intersectEdgesWithinFace(face) {
        const intersections = [];
        const edges = face.getEdges();
        for (let i = 0; i < edges.length; i++) {
            for (let j = i + 1; j < edges.length; j++) {
                const edgeIntersections = this.intersectEdges(edges[i], edges[j]);
                if (edgeIntersections.length > 0) {
                    intersections.push(...edgeIntersections);
                }
            }
        }
        return intersections;
    }
    findEdgeParameter(edge, point) {
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const curvePoint = edge.evaluate(t);
            if (VectorUtils.isClose(curvePoint, point)) {
                return t;
            }
        }
        return null;
    }
    findSurfaceParameterAtPoint(surface, point) {
        const tolerance = 0.0001;
        const stepU = 0.1;
        const stepV = 0.1;
        for (let u = 0; u <= 1; u += stepU) {
            for (let v = 0; v <= 1; v += stepV) {
                const pointOnSurface = surface.evaluate(u, v);
                if (VectorUtils.isClose(point, pointOnSurface, tolerance)) {
                    return {
                        u,
                        v
                    };
                }
            }
        }
        return null;
    }
    handleFaceIntersections() {
        const intersectionInfos = [];
        for (let i = 0; i < this.faces.length; i++) {
            for (let j = i + 1; j < this.faces.length; j++) {
                const faceIntersections = this.intersectFaces(this.faces[i], this.faces[j]);
                if (faceIntersections.length > 0) {
                    intersectionInfos.push({
                        face1: this.faces[i],
                        face2: this.faces[j],
                        intersections: faceIntersections
                    });
                }
            }
        }
        return intersectionInfos;
    }
    closestIntersection(otherBREP) {
        const closestIntersections = [];
        this.faces.forEach(faceA => {
            otherBREP.faces.forEach(faceB => {
                const intersections = this.intersectFaces(faceA, faceB);
                intersections.forEach(({points}) => {
                    points.forEach(point => {
                        const distance = VectorUtils.distanceBetweenPoints(point, faceA.evaluate(0.5));
                        closestIntersections.push({
                            point,
                            distance
                        });
                    });
                });
            });
        });
        if (closestIntersections.length === 0)
            return null;
        return closestIntersections.reduce((min, current) => current.distance < min.distance ? current : min);
    }
    intersectsWith(otherBREP) {
        const intersections = [];
        this.faces.forEach(faceA => {
            otherBREP.faces.forEach(faceB => {
                const faceIntersections = this.intersectFaces(faceA, faceB);
                if (faceIntersections.length > 0) {
                    intersections.push({
                        faceA,
                        faceB,
                        intersections: faceIntersections
                    });
                }
            });
        });
        return intersections;
    }
}
class Face {
    constructor(nurbsSurface) {
        this.nurbsSurface = nurbsSurface;
        this.edges = [];
    }
    addEdge(edge) {
        this.edges.push(edge);
    }
    getEdges() {
        return this.edges;
    }
    getArea() {
        const vertices = this.getVertices();
        let area = 0;
        for (let i = 1; i < vertices.length - 1; i++) {
            area += this.triangleArea(vertices[0], vertices[i], vertices[i + 1]);
        }
        return area;
    }
    getNormal() {
        const vertices = this.getVertices();
        if (vertices.length < 3) {
            throw new Error('Not enough vertices to calculate the normal.');
        }
        const v1 = VectorUtils.subtractVectors(vertices[1].getCoordinates(), vertices[0].getCoordinates());
        const v2 = VectorUtils.subtractVectors(vertices[2].getCoordinates(), vertices[0].getCoordinates());
        return VectorUtils.normalizeVector(VectorUtils.crossProduct(v1, v2));
    }
    getEdgePoints() {
        return this.edges.map(edge => edge.nurbsCurve.controlPoints);
    }
    removeEdge(edge) {
        const index = this.edges.indexOf(edge);
        if (index !== -1) {
            this.edges.splice(index, 1);
        } else {
            throw new Error('Edge does not exist in the face.');
        }
    }
    getVertices() {
        const verticesSet = new Set();
        this.edges.forEach(edge => {
            edge.getVertices().forEach(controlPoint => {
                verticesSet.add(controlPoint.toString());
            });
        });
        return Array.from(verticesSet).map(vertexString => vertexString.split(',').map(Number));
    }
    isClose(vertex1, vertex2, tolerance = 0.000001) {
        return vertex1.every((val, index) => Math.abs(val - vertex2[index]) < tolerance);
    }
    calculateArea() {
        // Calculates the face area using triangulation
        const vertices = this.getVertices();
        let area = 0;
        for (let i = 1; i < vertices.length - 1; i++) {
            area += this.triangleArea(vertices[0], vertices[i], vertices[i + 1]);
        }
        return area;
    }
    calculateCentroid() {
        // Computes the centroid of the face by averaging the vertices
        const vertices = this.getVertices();
        const total = vertices.reduce((acc, vertex) => VectorUtils.addVectors(acc, vertex), [
            0,
            0,
            0
        ]);
        return total.map(coord => coord / vertices.length);
    }
    isPointInside(point) {
        const edges = this.getEdges();
        let intersections = 0;
        for (const edge of edges) {
            const edgePoints = edge.nurbsCurve.controlPoints;
            if (this.isClose(edgePoints[0], point) || this.isClose(edgePoints[1], point)) {
                return true;
            }
            if (this.rayIntersectsEdge(edge, point)) {
                intersections++;
            }
        }
        return intersections % 2 === 1;
    }
    calculateNormal() {
        const vertices = this.getVertices();
        if (vertices.length < 3) {
            throw new Error('Not enough vertices to calculate the normal.');
        }
        const v1 = VectorUtils.subtractVectors(vertices[1].getCoordinates(), vertices[0].getCoordinates());
        const v2 = VectorUtils.subtractVectors(vertices[2].getCoordinates(), vertices[0].getCoordinates());
        return VectorUtils.normalizeVector(VectorUtils.crossProduct(v1, v2));
    }
    removeEdgeByControlPoints(controlPoint1, controlPoint2) {
        const edge = this.edges.find(e => this.isClose(e.nurbsCurve.controlPoints[0], controlPoint1) && this.isClose(e.nurbsCurve.controlPoints[1], controlPoint2));
        if (edge) {
            this.removeEdge(edge);
        } else {
            throw new Error('Edge not found for removal.');
        }
    }
    rayIntersectsEdge(edge, point) {
        const [p1, p2] = edge.nurbsCurve.controlPoints;
        const t = (point[1] - p1[1]) / (p2[1] - p1[1]);
        const xIntersection = p1[0] + t * (p2[0] - p1[0]);
        return xIntersection >= point[0];
    }
    triangleArea(p0, p1, p2) {
        const a = this.distanceBetweenPoints(p0, p1);
        const b = this.distanceBetweenPoints(p1, p2);
        const c = this.distanceBetweenPoints(p2, p0);
        const s = (a + b + c) / 2;
        return Math.sqrt(s * (s - a) * (s - b) * (s - c));
    }
    clearEdges() {
        this.edges = [];
    }
    isFaceValid() {
        return this.getEdges().length > 0 && this.getEdges().every(edge => edge.isValid());
    }
    isPlanar(tolerance = 0.000001) {
        const normalVector = this.getNormal();
        const referencePoint = this.edges[0].nurbsCurve.controlPoints[0];
        return this.edges.every(edge => {
            return edge.getVertices().every(vertex => {
                const vector = VectorUtils.subtractVectors(vertex.getCoordinates(), referencePoint);
                const distance = Math.abs(normalVector.reduce((sum, val, index) => val * vector[index] + sum, 0));
                return distance < tolerance;
            });
        });
    }
    distanceBetweenPoints(p1, p2) {
        return Math.sqrt(p1.reduce((sum, val, index) => sum + (val - p2[index]) ** 2, 0));
    }
    intersectEdgesWithinFace() {
        const intersections = [];
        const edges = this.getEdges();
        for (let i = 0; i < edges.length; i++) {
            for (let j = i + 1; j < edges.length; j++) {
                const intersectionPoints = edges[i].intersectEdges(edges[j]);
                if (intersectionPoints.length > 0) {
                    intersections.push(...intersectionPoints);
                }
            }
        }
        return intersections;
    }
    intersectEdgeWithPlane(edge, planePoint, planeNormal) {
        const controlPoints = edge.getControlPoints();
        for (let i = 0; i < controlPoints.length - 1; i++) {
            const p0 = controlPoints[i];
            const p1 = controlPoints[i + 1];
            const direction = VectorUtils.subtractVectors(p1, p0);
            const dotProductDirectionNormal = planeNormal.reduce((sum, val, index) => sum + val * direction[index], 0);
            if (Math.abs(dotProductDirectionNormal) > Number.EPSILON) {
                const t = planeNormal.reduce((sum, val, index) => sum + val * (planePoint[index] - p0[index]), 0) / dotProductDirectionNormal;
                if (t >= 0 && t <= 1) {
                    return p0.map((coord, index) => coord + t * direction[index]);
                }
            }
        }
        return null;
    }
    intersectWith(face) {
        const intersections = [];
        this.edges.forEach(edge1 => {
            face.edges.forEach(edge2 => {
                const edgeIntersections = edge1.intersectEdges(edge2);
                if (edgeIntersections.length > 0) {
                    intersections.push({
                        edge1,
                        edge2,
                        points: edgeIntersections
                    });
                }
            });
        });
        return intersections;
    }
    intersectEdgesWithOtherFace(otherFace) {
        const intersections = [];
        this.getEdges().forEach(edge1 => {
            otherFace.getEdges().forEach(edge2 => {
                const edgeIntersections = edge1.intersectEdges(edge2);
                if (edgeIntersections.length > 0) {
                    intersections.push({
                        edge1,
                        edge2,
                        points: edgeIntersections
                    });
                }
            });
        });
        return intersections;
    }
}
class Edge {
    constructor(nurbsCurve) {
        this.nurbsCurve = nurbsCurve;
    }
    evaluate(t) {
        return this.nurbsCurve.evaluate(t);
    }
    length() {
        // Calculates the cumulative length of the edge
        let totalLength = 0;
        const points = this.nurbsCurve.getPoints(0.01);
        for (let i = 0; i < points.length - 1; i++) {
            totalLength += Edge.distanceBetweenPoints(points[i], points[i + 1]);
        }
        return totalLength;
    }
    getMidPoint() {
        // Returns the midpoint of the edge curve by evaluating at t=0.5
        return this.nurbsCurve.evaluate(0.5);
    }
    getControlPoints() {
        // Retrieves the control points of the edge's NURBS curve
        return this.nurbsCurve.controlPoints;
    }
    splitAtParameter(t) {
        const controlPoints = this.nurbsCurve.controlPoints;
        const splitPoint = this.evaluate(t);
        const leftPoints = controlPoints.slice(0, Math.floor(controlPoints.length / 2)).concat(splitPoint);
        const rightPoints = [splitPoint].concat(controlPoints.slice(Math.floor(controlPoints.length / 2)));
        const leftCurve = new NURBS_Curve(this.nurbsCurve.degree, leftPoints, this.nurbsCurve.knots);
        const rightCurve = new NURBS_Curve(this.nurbsCurve.degree, rightPoints, this.nurbsCurve.knots);
        return {
            leftEdge: new Edge(leftCurve),
            rightEdge: new Edge(rightCurve)
        };
    }
    getVertices() {
        return this.nurbsCurve.controlPoints.map(cp => new Vertex(cp));
    }
    getLength() {
        return this.nurbsCurve.length();
    }
    setControlPoints(newControlPoints) {
        this.nurbsCurve.controlPoints = newControlPoints;
    }
    clearControlPoints() {
        this.nurbsCurve.controlPoints = [];
    }
    isValid() {
        return this.nurbsCurve !== null && this.nurbsCurve.controlPoints.length > 1;
    }
    isStraightLine(tolerance = 0.000001) {
        const [startPoint, endPoint] = [
            this.nurbsCurve.controlPoints[0],
            this.nurbsCurve.controlPoints[this.nurbsCurve.controlPoints.length - 1]
        ];
        const lineVector = VectorUtils.subtractVectors(endPoint, startPoint);
        const normalizedLine = VectorUtils.normalizeVector(lineVector);
        return this.nurbsCurve.controlPoints.every(point => {
            const vector = VectorUtils.subtractVectors(point, startPoint);
            const distance = Math.abs(normalizedLine.reduce((sum, val, index) => val * vector[index] + sum, 0));
            return distance < tolerance;
        });
    }
    static distanceBetweenPoints(p1, p2) {
        return Math.sqrt(p1.reduce((sum, val, idx) => sum + (val - p2[idx]) ** 2, 0));
    }
    intersectEdges(otherEdge) {
        const intersections = [];
        const stepSize = 0.01;
        for (let t1 = 0; t1 <= 1; t1 += stepSize) {
            const point1 = this.evaluate(t1);
            for (let t2 = 0; t2 <= 1; t2 += stepSize) {
                const point2 = otherEdge.evaluate(t2);
                if (VectorUtils.isClose(point1, point2)) {
                    intersections.push({
                        point: point1,
                        param1: t1,
                        param2: t2
                    });
                }
            }
        }
        return intersections;
    }
    findParameterAtPoint(point) {
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const curvePoint = this.evaluate(t);
            if (VectorUtils.isClose(curvePoint, point)) {
                return t;
            }
        }
        return null;
    }
    intersectWithOtherEdge(otherEdge) {
        return BREP.prototype.intersectEdges(this, otherEdge);
    }
    splitAtIntersections(otherEdge) {
        const intersections = this.intersectWithOtherEdge(otherEdge);
        if (!intersections || intersections.length === 0) {
            return [this];
        }
        const splitEdges = [];
        let remainingCurve = this.nurbsCurve;
        intersections.forEach(({param1}) => {
            const splitResult = remainingCurve.splitAtParameter(param1);
            splitEdges.push(new Edge(splitResult.leftCurve));
            remainingCurve = splitResult.rightCurve;
        });
        splitEdges.push(new Edge(remainingCurve));
        return splitEdges;
    }
    intersectWith(otherEdge) {
        return BREP.prototype.intersectEdges(this, otherEdge);
    }
}
class Vertex {
    constructor(point) {
        this.point = point;
    }
    getCoordinates() {
        return this.point;
    }
    isEqualTo(otherVertex) {
        // Checks vertex equality by comparing coordinates
        return this.getCoordinates().every((coord, index) => coord === otherVertex.getCoordinates()[index]);
    }
    calculateNormal() {
        const faces = this.getFaces();
        const normals = faces.map(face => face.calculateNormal());
        const average = normals.reduce((sum, normal) => VectorUtils.addVectors(sum, normal), [
            0,
            0,
            0
        ]);
        return VectorUtils.normalizeVector(average);
    }
    getFaces() {
        const facesConnected = [];
        const edges = BREP.getEdgesByVertex(this);
        edges.forEach(edge => {
            facesConnected.push(BREP.findFaceByEdge(edge));
        });
        return facesConnected;
    }
    clear() {
        this.point = null;
    }
}
class BREP_Kernel {
    constructor() {
        this.BREP = new BREP();
    }
    addFace(nurbsSurface) {
        return this.BREP.addFace(nurbsSurface);
    }
    addEdge(face, nurbsCurve) {
        const edge = new Edge(nurbsCurve);
        face.addEdge(edge);
        return edge;
    }
    getVerticesForFace(face) {
        // Retrieves vertices for a specified face without duplicates
        const vertices = [];
        for (const edge of face.getEdges()) {
            const controlPoints = edge.nurbsCurve.controlPoints;
            for (const point of controlPoints) {
                const vertex = new Vertex(point);
                if (!vertices.some(v => this.isClose(v.getCoordinates(), vertex.getCoordinates()))) {
                    vertices.push(vertex);
                }
            }
        }
        return vertices;
    }
    isClose(point1, point2, tolerance = 0.000001) {
        return point1.every((val, index) => Math.abs(val - point2[index]) < tolerance);
    }
    removeFace(face) {
        this.BREP.removeFace(face);
    }
    removeEdge(face, edge) {
        const edges = face.getEdges();
        const index = edges.indexOf(edge);
        if (index !== -1) {
            edges.splice(index, 1);
        } else {
            throw new Error('Edge does not exist in the face.');
        }
    }
    updateFace(face, newNurbsSurface) {
        face.nurbsSurface = newNurbsSurface;
    }
    updateEdge(edge, newNurbsCurve) {
        const existingEdge = this.BREP.getEdges().find(e => e === edge);
        if (existingEdge) {
            existingEdge.nurbsCurve = newNurbsCurve;
        } else {
            throw new Error('Edge not found for update.');
        }
    }
    mergeEdges(edge1, edge2) {
        if (!this.edgesShareVertices(edge1, edge2)) {
            throw new Error('Cannot merge edges that do not share vertices.');
        }
        const combinedControlPoints = [
            ...edge1.nurbsCurve.controlPoints.slice(0, -1),
            ...edge2.nurbsCurve.controlPoints
        ];
        const degree = Math.max(edge1.nurbsCurve.degree, edge2.nurbsCurve.degree);
        const newKnots = this.mergeKnotVectors(edge1.nurbsCurve.knots, edge2.nurbsCurve.knots);
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    combineCurves(curve1, curve2) {
        const combinedControlPoints = [
            ...curve1.controlPoints,
            ...curve2.controlPoints.slice(1)
        ];
        const degree = Math.max(curve1.degree, curve2.degree);
        const newKnots = this.mergeKnotVectors(curve1.knots.concat(curve2.knots));
        return new NURBS_Curve(degree, combinedControlPoints, newKnots);
    }
    mergeKnotVectors(knots) {
        return Array.from(new Set(knots.flat())).sort((a, b) => a - b);
    }
    cleanUpBREP() {
        for (const face of this.BREP.getFaces()) {
            face.edges = face.getEdges().filter(edge => edge.isValid());
        }
    }
    getEdgeVertices(edge) {
        return edge.nurbsCurve.controlPoints.map(cp => new Vertex(cp));
    }
    findSharedEdges(face1, face2) {
        const edgesFace1 = face1.getEdges();
        const edgesFace2 = face2.getEdges();
        return edgesFace1.filter(edge1 => edgesFace2.some(edge2 => edge1.nurbsCurve === edge2.nurbsCurve));
    }
    mergeFaces(face1, face2) {
        const sharedEdges = this.findSharedEdges(face1, face2);
        if (sharedEdges.length === 0) {
            throw new Error('No shared edges between the faces, cannot merge.');
        }
        const combinedControlPoints = [
            ...face1.nurbsSurface.controlPoints,
            ...face2.nurbsSurface.controlPoints
        ];
        const degreeU = Math.max(face1.nurbsSurface.degreeU, face2.nurbsSurface.degreeU);
        const degreeV = Math.max(face1.nurbsSurface.degreeV, face2.nurbsSurface.degreeV);
        const newKnotsU = this.mergeKnotVectors(face1.nurbsSurface.knotsU, face2.nurbsSurface.knotsU);
        const newKnotsV = this.mergeKnotVectors(face1.nurbsSurface.knotsV, face2.nurbsSurface.knotsV);
        const mergedSurface = new NURBS_Surface(degreeU, degreeV, combinedControlPoints, newKnotsU, newKnotsV);
        return new Face(mergedSurface);
    }
    splitFace(face, curve) {
        const edges = face.getEdges();
        const splitEdges = edges.filter(edge => curve.isClose(edge.nurbsCurve));
        const newControlPoints1 = [];
        const newControlPoints2 = [];
        splitEdges.forEach(edge => {
            const controlPoints = edge.nurbsCurve.controlPoints;
            const midPointIndex = Math.floor(controlPoints.length / 2);
            newControlPoints1.push(...controlPoints.slice(0, midPointIndex));
            newControlPoints2.push(...controlPoints.slice(midPointIndex));
        });
        const newSurface1 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints1, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newSurface2 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints2, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newFace1 = new Face(newSurface1);
        const newFace2 = new Face(newSurface2);
        this.BREP.removeFace(face);
        this.BREP.addFace(newFace1);
        this.BREP.addFace(newFace2);
        return {
            newFace1,
            newFace2
        };
    }
    getSharedVertices(face1, face2) {
        const verticesFace1 = this.getVerticesForFace(face1);
        const verticesFace2 = this.getVerticesForFace(face2);
        return verticesFace1.filter(v1 => verticesFace2.some(v2 => this.isClose(v1.getCoordinates(), v2.getCoordinates())));
    }
    edgesShareVertices(edge1, edge2) {
        const verticesEdge1 = this.getEdgeVertices(edge1);
        const verticesEdge2 = this.getEdgeVertices(edge2);
        return verticesEdge1.some(v1 => verticesEdge2.some(v2 => this.isClose(v1.getCoordinates(), v2.getCoordinates())));
    }
    scaleBREP(factor) {
        this.BREP.getFaces().forEach(face => {
            face.getEdges().forEach(edge => {
                edge.nurbsCurve.controlPoints = edge.nurbsCurve.controlPoints.map(cp => cp.map(coord => coord * factor));
            });
        });
    }
    translateBREP(vector) {
        this.BREP.getFaces().forEach(face => {
            face.getEdges().forEach(edge => {
                edge.nurbsCurve.controlPoints = edge.nurbsCurve.controlPoints.map(cp => cp.map((coord, index) => coord + vector[index]));
            });
        });
    }
    addVertex(face, point) {
        const vertex = new Vertex(point);
        if (!face.getVertices().some(v => this.isClose(v.getCoordinates(), vertex.getCoordinates()))) {
            face.addEdge(vertex);
        }
        return vertex;
    }
    removeVertex(vertex) {
        this.BREP.removeVertexFromEdges(vertex);
    }
    getVertexByCoordinates(face, coordinates) {
        const vertices = face.getVertices();
        return vertices.find(vertex => this.isClose(vertex.getCoordinates(), coordinates)) || null;
    }
    isEdgeShared(edge1, edge2) {
        const verticesEdge1 = this.getEdgeVertices(edge1);
        const verticesEdge2 = this.getEdgeVertices(edge2);
        return verticesEdge1.some(v1 => verticesEdge2.some(v2 => this.isClose(v1.getCoordinates(), v2.getCoordinates())));
    }
    isFaceShared(face1, face2) {
        const sharedEdges = this.findSharedEdges(face1, face2);
        return sharedEdges.length > 0;
    }
    findSharedVertices(face1, face2) {
        const verticesFace1 = this.getVerticesForFace(face1);
        const verticesFace2 = this.getVerticesForFace(face2);
        return verticesFace1.filter(v1 => verticesFace2.some(v2 => this.isClose(v1.getCoordinates(), v2.getCoordinates())));
    }
    mergeSharedFaces() {
        const faces = this.BREP.getFaces();
        const mergedFaces = [];
        for (let i = 0; i < faces.length; i++) {
            for (let j = i + 1; j < faces.length; j++) {
                const sharedEdges = this.findSharedEdges(faces[i], faces[j]);
                if (sharedEdges.length > 0) {
                    const mergedFace = this.mergeFaces(faces[i], faces[j]);
                    mergedFaces.push(mergedFace);
                    this.removeFace(faces[i]);
                    this.removeFace(faces[j]);
                    break;
                }
            }
        }
        mergedFaces.forEach(face => this.BREP.addFace(face));
    }
    splitFaceByPlane(face, point, normal) {
        const edges = face.getEdges();
        const aboveEdges = [];
        const belowEdges = [];
        edges.forEach(edge => {
            const edgeMidPoint = edge.getMidPoint();
            if (this.isPointAbovePlane(edgeMidPoint, point, normal)) {
                aboveEdges.push(edge);
            } else {
                belowEdges.push(edge);
            }
        });
        const aboveFace = new Face(face.nurbsSurface);
        aboveEdges.forEach(edge => aboveFace.addEdge(edge));
        const belowFace = new Face(face.nurbsSurface);
        belowEdges.forEach(edge => belowFace.addEdge(edge));
        this.removeFace(face);
        this.BREP.addFace(aboveFace);
        this.BREP.addFace(belowFace);
    }
    isPointAbovePlane(point, planePoint, normal) {
        const vectorToPoint = VectorUtils.subtractVectors(point, planePoint);
        const dotProduct = vectorToPoint.reduce((sum, val, index) => sum + val * normal[index], 0);
        return dotProduct > 0;
    }
    findClosestVertex(point) {
        let closestVertex = null;
        let minDistance = Infinity;
        for (const face of this.BREP.getFaces()) {
            const vertices = this.getVerticesForFace(face);
            for (const vertex of vertices) {
                const distance = this.distanceBetweenPoints(vertex.getCoordinates(), point);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestVertex = vertex;
                }
            }
        }
        return closestVertex;
    }
    distanceBetweenPoints(p1, p2) {
        return Math.sqrt(p1.reduce((sum, val, idx) => sum + (val - p2[idx]) ** 2, 0));
    }
    reverseEdges(face) {
        const edges = face.getEdges();
        for (const edge of edges) {
            edge.nurbsCurve.controlPoints.reverse();
            edge.nurbsCurve.knots.reverse();
        }
    }
    getFaceByEdge(edge) {
        for (const face of this.BREP.getFaces()) {
            if (face.getEdges().includes(edge)) {
                return face;
            }
        }
        return null;
    }
    getEdgeByVertices(vertex1, vertex2) {
        for (const face of this.BREP.getFaces()) {
            const edges = face.getEdges();
            for (const edge of edges) {
                if (this.isClose(edge.nurbsCurve.controlPoints[0], vertex1.getCoordinates()) && this.isClose(edge.nurbsCurve.controlPoints[1], vertex2.getCoordinates()) || this.isClose(edge.nurbsCurve.controlPoints[0], vertex2.getCoordinates()) && this.isClose(edge.nurbsCurve.controlPoints[1], vertex1.getCoordinates())) {
                    return edge;
                }
            }
        }
        return null;
    }
    splitEdge(edge, parameter) {
        const controlPoints = edge.nurbsCurve.controlPoints;
        const splitPoint = edge.evaluate(parameter);
        const newControlPoints = [...controlPoints];
        newControlPoints.splice(1, 0, splitPoint);
        const newCurve = new NURBS_Curve(edge.nurbsCurve.degree, newControlPoints, edge.nurbsCurve.knots);
        return newCurve;
    }
    getFaceVertices(face) {
        return face.getVertices();
    }
    updateEdgeControlPoints(edge, newControlPoints) {
        edge.nurbsCurve.controlPoints = newControlPoints;
    }
    findFaceByVertex(vertex) {
        for (const face of this.BREP.getFaces()) {
            const vertices = face.getVertices();
            if (vertices.some(v => this.isClose(v.getCoordinates(), vertex.getCoordinates()))) {
                return face;
            }
        }
        return null;
    }
    isEdgeOnFace(edge, face) {
        return face.getEdges().includes(edge);
    }
    mergeFacesByVertices(vertex1, vertex2) {
        const face1 = this.getFaceByVertices(vertex1, vertex2);
        const face2 = this.getFaceByVertices(vertex2, vertex1);
        if (!face1 || !face2) {
            throw new Error('No faces are sharing the specified vertices.');
        }
        return this.mergeFaces(face1, face2);
    }
    splitEdgeByVertices(vertex1, vertex2) {
        const edge = this.findEdgeWithVertices(vertex1, vertex2);
        if (!edge) {
            throw new Error('No edge found between the specified vertices.');
        }
        const midPoint = edge.getMidPoint();
        const newCurve = this.splitEdge(edge, 0.5);
        // Split at midPoint
        return newCurve;
    }
    splitFaceByCurve(face, curve) {
        const edges = face.getEdges();
        const splitEdges = edges.filter(edge => curve.isClose(edge.nurbsCurve));
        const newControlPoints1 = [];
        const newControlPoints2 = [];
        splitEdges.forEach(edge => {
            const controlPoints = edge.nurbsCurve.controlPoints;
            const midPointIndex = Math.floor(controlPoints.length / 2);
            newControlPoints1.push(...controlPoints.slice(0, midPointIndex));
            newControlPoints2.push(...controlPoints.slice(midPointIndex));
        });
        const newSurface1 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints1, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newSurface2 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints2, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newFace1 = new Face(newSurface1);
        const newFace2 = new Face(newSurface2);
        this.removeFace(face);
        this.addFace(newFace1);
        this.addFace(newFace2);
        return [
            newFace1,
            newFace2
        ];
    }
    updateVertex(vertex, newCoordinates) {
        const foundVertex = this.findVertex(vertex);
        if (foundVertex) {
            foundVertex.point = newCoordinates;
            this.updateConnectedEdges(foundVertex);
        } else {
            throw new Error('Vertex not found.');
        }
    }
    findEdgeWithVertices(vertex1, vertex2) {
        const edges = this.BREP.getEdges();
        return edges.find(edge => {
            const controlPoints = edge.nurbsCurve.controlPoints;
            return this.isClose(controlPoints[0], vertex1.getCoordinates()) && this.isClose(controlPoints[1], vertex2.getCoordinates()) || this.isClose(controlPoints[0], vertex2.getCoordinates()) && this.isClose(controlPoints[1], vertex1.getCoordinates());
        }) || null;
    }
    findEdgeByControlPoints(controlPoint1, controlPoint2) {
        const edges = this.BREP.getEdges();
        return edges.find(edge => {
            const controlPoints = edge.getControlPoints();
            return this.isClose(controlPoints[0], controlPoint1) && this.isClose(controlPoints[1], controlPoint2) || this.isClose(controlPoints[0], controlPoint2) && this.isClose(controlPoints[1], controlPoint1);
        }) || null;
    }
    isFaceValid(face) {
        return face.getEdges().length > 0 && face.getEdges().every(edge => edge.isValid());
    }
    splitFaceByEdge(face, edge) {
        const edges = face.getEdges();
        const splitEdges = edges.filter(e => e === edge);
        const newControlPoints1 = [];
        const newControlPoints2 = [];
        splitEdges.forEach(edge => {
            const controlPoints = edge.nurbsCurve.controlPoints;
            const midPointIndex = Math.floor(controlPoints.length / 2);
            newControlPoints1.push(...controlPoints.slice(0, midPointIndex + 1));
            newControlPoints2.push(...controlPoints.slice(midPointIndex));
        });
        const newSurface1 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints1, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newSurface2 = new NURBS_Surface(face.nurbsSurface.degreeU, face.nurbsSurface.degreeV, newControlPoints2, face.nurbsSurface.knotsU, face.nurbsSurface.knotsV);
        const newFace1 = new Face(newSurface1);
        const newFace2 = new Face(newSurface2);
        this.removeFace(face);
        this.addFace(newFace1);
        this.addFace(newFace2);
        return [
            newFace1,
            newFace2
        ];
    }
    calculateFaceArea(face) {
        return face.getEdges().reduce((totalArea, edge) => totalArea + edge.length(), 0);
    }
    combineEdges(edgeList) {
        if (!edgeList || edgeList.length < 2) {
            throw new Error('At least two edges are required for combination.');
        }
        const combinedControlPoints = edgeList.flatMap(edge => edge.nurbsCurve.controlPoints);
        const degree = Math.max(...edgeList.map(edge => edge.nurbsCurve.degree));
        const newKnots = this.mergeKnotVectors(edgeList.map(edge => edge.nurbsCurve.knots));
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    mergeFacesByEdge(face1, face2) {
        const sharedEdges = this.findSharedEdges(face1, face2);
        if (sharedEdges.length === 0) {
            throw new Error('No shared edges between the faces, cannot merge.');
        }
        const combinedControlPoints = sharedEdges.reduce((points, edge) => [
            ...points,
            ...edge.nurbsCurve.controlPoints
        ], []);
        const degree = Math.max(face1.nurbsSurface.degreeU, face2.nurbsSurface.degreeU);
        const newKnotsU = this.mergeKnotVectors(face1.nurbsSurface.knotsU, face2.nurbsSurface.knotsU);
        const newKnotsV = this.mergeKnotVectors(face1.nurbsSurface.knotsV, face2.nurbsSurface.knotsV);
        const mergedSurface = new NURBS_Surface(degree, degree, combinedControlPoints, newKnotsU, newKnotsV);
        return new Face(mergedSurface);
    }
    removeEdgeFromFace(face, edge) {
        const edges = face.getEdges();
        const index = edges.indexOf(edge);
        if (index !== -1) {
            edges.splice(index, 1);
        } else {
            throw new Error('Edge does not exist in the face.');
        }
    }
    mergeEdgesByFace(face1, face2) {
        const edgesFace1 = this.getEdgesByFace(face1);
        const edgesFace2 = this.getEdgesByFace(face2);
        const sharedEdges = edgesFace1.filter(edge1 => edgesFace2.some(edge2 => this.isEdgeShared(edge1, edge2)));
        if (sharedEdges.length === 0) {
            throw new Error('No shared edges between the faces, cannot merge.');
        }
        const combinedControlPoints = sharedEdges.reduce((points, edge) => {
            return [
                ...points,
                ...edge.nurbsCurve.controlPoints
            ];
        }, []);
        const degree = Math.max(...sharedEdges.map(edge => edge.nurbsCurve.degree));
        const newKnots = this.mergeKnotVectors(sharedEdges.map(edge => edge.nurbsCurve.knots));
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    isFaceConnected(face1, face2) {
        const edgesFace1 = this.getEdgesByFace(face1);
        const edgesFace2 = this.getEdgesByFace(face2);
        return edgesFace1.some(edge1 => edgesFace2.some(edge2 => this.isEdgeShared(edge1, edge2)));
    }
    getAdjacentFaces(face) {
        return this.BREP.getFaces().filter(f => f.getEdges().some(edge => face.getEdges().includes(edge)));
    }
    findFacesByEdge(edge) {
        return this.BREP.getFaces().filter(face => face.getEdges().includes(edge));
    }
    splitFaceByPoint(face, point) {
        const edges = face.getEdges();
        const splitEdges = [];
        edges.forEach(edge => {
            const parameter = this.findParameterAtPoint(edge.nurbsCurve, point);
            if (parameter !== null) {
                splitEdges.push(edge.splitAtParameter(parameter));
            }
        });
        const newFace = new Face(face.nurbsSurface);
        splitEdges.forEach(edge => newFace.addEdge(edge));
        this.removeFace(face);
        this.addFace(newFace);
        return newFace;
    }
    findParameterAtPoint(curve, point) {
        const stepSize = 0.01;
        for (let t = 0; t <= 1; t += stepSize) {
            const curvePoint = curve.evaluate(t);
            if (this.isClose(curvePoint, point)) {
                return t;
            }
        }
        return null;
    }
    getFaceArea(face) {
        return face.getArea();
    }
    addEdgeToFace(face, nurbsCurve) {
        const edge = new Edge(nurbsCurve);
        face.addEdge(edge);
        return edge;
    }
    calculateSurfaceArea() {
        return this.BREP.getFaces().reduce((totalArea, face) => totalArea + face.getArea(), 0);
    }
    findVertexByCoordinates(face, coordinates) {
        const vertices = face.getVertices();
        return vertices.find(vertex => this.isClose(vertex.getCoordinates(), coordinates)) || null;
    }
    addVertexToFace(face, point) {
        const vertex = new Vertex(point);
        if (!face.getVertices().some(v => this.isClose(v.getCoordinates(), vertex.getCoordinates()))) {
            face.addEdge(vertex);
        }
        return vertex;
    }
    removeVertexFromFace(face, vertex) {
        face.removeEdgeByControlPoints(vertex.getCoordinates(), vertex.getCoordinates());
    }
    weldEdges(edge1, edge2) {
        if (!this.edgesShareVertices(edge1, edge2)) {
            throw new Error('Cannot weld edges that do not share vertices.');
        }
        const combinedControlPoints = [
            ...edge1.nurbsCurve.controlPoints.slice(0, -1),
            ...edge2.nurbsCurve.controlPoints
        ];
        const degree = Math.max(edge1.nurbsCurve.degree, edge2.nurbsCurve.degree);
        const newKnots = this.mergeKnotVectors(edge1.nurbsCurve.knots, edge2.nurbsCurve.knots);
        return new Edge(new NURBS_Curve(degree, combinedControlPoints, newKnots));
    }
    removeOrphanedEdges() {
        const orphanedEdges = this.BREP.getEdges().filter(edge => {
            return !this.BREP.getFaces().some(face => face.getEdges().includes(edge));
        });
        orphanedEdges.forEach(edge => {
            this.BREP.getFaces().forEach(face => {
                face.removeEdge(edge);
            });
        });
    }
    validate() {
        const faceValidity = this.BREP.getFaces().every(face => face.isFaceValid());
        const edgeValidity = this.BREP.getEdges().every(edge => edge.isValid());
        return faceValidity && edgeValidity;
    }
    computeFaceNormals() {
        this.BREP.getFaces().forEach(face => {
            const edges = face.getEdges();
            if (edges.length < 3) {
                throw new Error('A face needs at least three edges to define a normal.');
            }
            const normal = face.calculateNormal();
            face.normal = normal;
        });
    }
    clearBREP() {
        this.BREP.faces.forEach(face => face.clearEdges());
        this.BREP.clearFaces();
    }
    findVertex(vertex) {
        return this.BREP.getVertices().find(v => this.isClose(v.getCoordinates(), vertex.getCoordinates())) || null;
    }
    isEdgeValid(edge) {
        return edge && edge.nurbsCurve.controlPoints.length > 1;
    }
    mergeEdgesByVertices(vertex1, vertex2) {
        const edge1 = this.findEdgeWithVertices(vertex1, vertex2);
        if (!edge1) {
            throw new Error('No edge found between the specified vertices.');
        }
        const combinedCurve = this.combineCurves(edge1.nurbsCurve, edge1.nurbsCurve);
        return new Edge(combinedCurve);
    }
    mergeCurves(curve1, curve2) {
        if (!this.edgesShareVertices(curve1, curve2)) {
            throw new Error('Cannot merge curves that do not share vertices.');
        }
        const combinedControlPoints = [
            ...curve1.controlPoints.slice(0, -1),
            ...curve2.controlPoints
        ];
        const degree = Math.max(curve1.degree, curve2.degree);
        const newKnots = this.mergeKnotVectors(curve1.knots, curve2.knots);
        return new NURBS_Curve(degree, combinedControlPoints, newKnots);
    }
    clear() {
        this.BREP.faces.forEach(face => face.clearEdges());
        this.BREP.clearFaces();
    }
    getEdgesByFace(face) {
        return face.getEdges();
    }
    findFaceByEdge(edge) {
        for (const face of this.BREP.getFaces()) {
            if (face.getEdges().includes(edge)) {
                return face;
            }
        }
        return null;
    }
    calculateVolume() {
        return this.BREP.getFaces().reduce((totalVolume, face) => totalVolume + face.getArea() * face.height, 0);
    }
    getEdges() {
        return this.BREP.getEdges();
    }
}