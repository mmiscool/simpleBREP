const defaultTolerance = 0.000001;
class Point {
    constructor(x, y, z, id) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.id = id;
    }
    distanceTo(point) {
        const dx = this.x - point.x;
        const dy = this.y - point.y;
        const dz = this.z - point.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
    transform(matrix) {
        const [x, y, z] = [
            this.x,
            this.y,
            this.z
        ];
        this.x = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3];
        this.y = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3];
        this.z = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3];
    }
    equals(point, tolerance = defaultTolerance) {
        return this.distanceTo(point) <= tolerance;
    }
    toArray() {
        return [
            this.x,
            this.y,
            this.z
        ];
    }
    lerp(point, t) {
        const x = this.x + t * (point.x - this.x);
        const y = this.y + t * (point.y - this.y);
        const z = this.z + t * (point.z - this.z);
        return new Point(x, y, z);
    }
    merge(vertex, tolerance) {
        const distance = this.distanceTo(vertex.point);
        return distance <= tolerance;
    }
    isCoincident(vertex, tolerance) {
        return this.distanceTo(vertex.point) <= tolerance;
    }
    clone() {
        return new Point(this.x, this.y, this.z, this.id);
    }
}
class Vertex {
    constructor(point, id) {
        this.point = point;
        this.id = id;
    }
    getPosition() {
        return this.point.toArray();
    }
    transform(matrix) {
        this.point.transform(matrix);
    }
    merge(vertex, tolerance) {
        return this.point.merge(vertex, tolerance);
    }
    isCoincident(vertex, tolerance) {
        return this.point.isCoincident(vertex.point, tolerance);
    }
    clone() {
        return new Vertex(this.point.clone(), this.id);
    }
    add(vertex) {
        const merged = this.merge(vertex, defaultTolerance);
        return merged ? this : new Vertex(vertex.point.clone(), vertex.id);
    }
}
class NURBSCurve {
    constructor(degree, controlPoints, weights, knots, isClosed) {
        this.degree = degree;
        this.controlPoints = controlPoints;
        this.weights = weights;
        this.knots = knots;
        this.isClosed = isClosed;
    }
    evaluate(u) {
        const n = this.controlPoints.length - 1;
        const p = this.degree;
        let result = new Point(0, 0, 0);
        const denominators = [];
        let denominator = 0;
        for (let i = 0; i <= n; i++) {
            const span = this.findSpan(n, p, u);
            const basis = this.basisFunction(span, u, p, this.knots);
            denominator += this.weights[i] * basis;
            const point = this.controlPoints[i];
            result.x += point.x * this.weights[i] * basis;
            result.y += point.y * this.weights[i] * basis;
            result.z += point.z * this.weights[i] * basis;
        }
        if (denominator !== 0) {
            result.x /= denominator;
            result.y /= denominator;
            result.z /= denominator;
        }
        return result;
    }
    tangent(u) {
        const n = this.controlPoints.length - 1;
        const p = this.degree;
        const span = this.findSpan(n, p, u);
        const basisDerivatives = this.basisFunctionDerivatives(span, u, p, this.knots);
        const tangent = new Point(0, 0, 0);
        const denominator = this.getDenominator(u, span);
        for (let i = 0; i <= n; i++) {
            tangent.x += basisDerivatives[1][i] * this.weights[i] * this.controlPoints[i].x;
            tangent.y += basisDerivatives[1][i] * this.weights[i] * this.controlPoints[i].y;
            tangent.z += basisDerivatives[1][i] * this.weights[i] * this.controlPoints[i].z;
        }
        if (denominator !== 0) {
            tangent.x /= denominator;
            tangent.y /= denominator;
            tangent.z /= denominator;
        }
        return tangent;
    }
    split(u) {
        const span = this.findSpan(this.controlPoints.length - 1, this.degree, u);
        const leftControlPoints = [];
        const rightControlPoints = [];
        const leftWeights = [];
        const rightWeights = [];
        const leftBasis = [];
        const rightBasis = [];
        const n = this.controlPoints.length - 1;
        // Calculate left and right control points and weights based on the knot span
        for (let i = 0; i <= span; i++) {
            const basis = this.basisFunction(span, u, this.degree, this.knots);
            leftControlPoints.push(this.controlPoints[i]);
            leftWeights.push(this.weights[i]);
            leftBasis.push(basis);
        }
        for (let i = span + 1; i <= n; i++) {
            const basis = this.basisFunction(i, u, this.degree, this.knots);
            rightControlPoints.push(this.controlPoints[i]);
            rightWeights.push(this.weights[i]);
            rightBasis.push(basis);
        }
        // Create new left and right curves
        const leftCurve = new NURBSCurve(this.degree, leftControlPoints, leftWeights, this.knots.slice(0, span + 1), this.isClosed);
        const rightCurve = new NURBSCurve(this.degree, rightControlPoints, rightWeights, this.knots.slice(span + 1), this.isClosed);
        return [
            leftCurve,
            rightCurve
        ];
    }
    transform(matrix) {
        const transformedControlPoints = this.controlPoints.map(point => {
            const [x, y, z] = [
                point.x,
                point.y,
                point.z
            ];
            return new Point(matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3], matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3], matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3]);
        });
        this.controlPoints = transformedControlPoints;
    }
    projectToPlane(plane) {
        const projectedControlPoints = this.controlPoints.map(point => {
            const projection = this.calculateProjection(point, plane);
            return new Point(projection.x, projection.y, projection.z);
        });
        this.controlPoints = projectedControlPoints;
    }
    refineKnots() {
        const newKnots = [];
        const {degree, knots} = this;
        const n = knots.length - degree - 1;
        for (let i = 0; i <= n; i++) {
            const knotSpan = knots[i + 1] - knots[i];
            newKnots.push(knots[i]);
            // Introducing additional knots based on current knot span
            if (i < n) {
                const additionalKnotsCount = Math.floor(Math.max(1, knotSpan * 10));
                // Add more knots for finer refinement
                for (let j = 1; j <= additionalKnotsCount; j++) {
                    const newKnot = knots[i] + j * (knotSpan / (additionalKnotsCount + 1));
                    newKnots.push(newKnot);
                }
            }
        }
        newKnots.push(knots[n + 1]);
        this.knots = newKnots;
    }
    getBoundingBox() {
        const minX = Math.min(...this.controlPoints.map(point => point.x));
        const minY = Math.min(...this.controlPoints.map(point => point.y));
        const minZ = Math.min(...this.controlPoints.map(point => point.z));
        const maxX = Math.max(...this.controlPoints.map(point => point.x));
        const maxY = Math.max(...this.controlPoints.map(point => point.y));
        const maxZ = Math.max(...this.controlPoints.map(point => point.z));
        return {
            min: new Point(minX, minY, minZ),
            max: new Point(maxX, maxY, maxZ)
        };
    }
    getLength(precision) {
        let length = 0;
        const step = 1 / Math.max(1, precision);
        const n = this.controlPoints.length - 1;
        for (let i = 0; i < n; i++) {
            const startU = this.knots[i + this.degree] + step;
            const endU = this.knots[i + this.degree + 1];
            for (let u = startU; u <= endU; u += step) {
                const point1 = this.evaluate(u - step);
                const point2 = this.evaluate(u);
                length += point1.distanceTo(point2);
            }
        }
        return length;
    }
    derive(u, order) {
        const n = this.controlPoints.length - 1;
        const p = this.degree;
        const span = this.findSpan(n, p, u);
        const basisDerivatives = this.basisFunctionDerivatives(span, u, p, this.knots);
        const derivatives = Array.from({ length: order + 1 }, () => new Point(0, 0, 0));
        const denominator = this.getDenominator(u, span);
        for (let j = 0; j <= order; j++) {
            for (let i = 0; i <= n; i++) {
                if (j === 0) {
                    derivatives[j].x += basisDerivatives[j][i] * this.weights[i] * this.controlPoints[i].x;
                    derivatives[j].y += basisDerivatives[j][i] * this.weights[i] * this.controlPoints[i].y;
                    derivatives[j].z += basisDerivatives[j][i] * this.weights[i] * this.controlPoints[i].z;
                } else {
                    derivatives[j].x += basisDerivatives[j][i] * this.weights[i] * this.controlPoints[i].x;
                    derivatives[j].y += basisDerivatives[j][i] * this.weights[i] * this.controlPoints[i].y;
                    derivatives[j].z += basisDerivatives[j][i] * this.weights[i] * this.controlPoints[i].z;
                }
            }
        }
        if (denominator !== 0) {
            for (let j = 0; j <= order; j++) {
                derivatives[j].x /= denominator;
                derivatives[j].y /= denominator;
                derivatives[j].z /= denominator;
            }
        }
        return derivatives;
    }
    getControlPolygon() {
        return this.controlPoints.map(point => point.clone());
    }
    findSpan(n, p, u) {
        if (u >= this.knots[n + 1]) {
            return n;
        }
        if (u <= this.knots[p]) {
            return 0;
        }
        let low = p;
        let high = n + 1;
        const mid = Math.floor((low + high) / 2);
        while (u < this.knots[mid] || u >= this.knots[mid + 1]) {
            if (u < this.knots[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    basisFunction(span, u, p, knotVector) {
        const leftNumerator = (u - knotVector[span]) / (knotVector[span + 1] - knotVector[span]);
        const rightNumerator = (knotVector[span + 1 + p] - u) / (knotVector[span + 1 + p] - knotVector[span + 1]);
        let left = 0;
        let right = 0;
        if (p === 0) {
            return u >= knotVector[span] && u < knotVector[span + 1] ? 1 : 0;
        } else {
            left = this.basisFunction(span - 1, u, p - 1, knotVector) * leftNumerator;
            right = this.basisFunction(span, u, p - 1, knotVector) * rightNumerator;
        }
        return left + right;
    }
    basisFunctionDerivatives(span, u, p, knotVector) {
        const derivatives = Array.from({ length: p + 1 }, () => 0);
        const left = Array(p + 1).fill(0);
        const right = Array(p + 1).fill(0);
        const temp = Array(p + 1).fill(0);
        derivatives[0] = 1;
        for (let j = 1; j <= p; j++) {
            // Compute left and right values for basis function derivatives
            left[j] = u - knotVector[span + 1 - j];
            right[j] = knotVector[span + j] - u;
            const saved = 0;
            for (let r = 0; r < j; r++) {
                const tempValue = derivatives[r] / (right[r + 1] + left[j - r]);
                derivatives[r] = saved + right[r + 1] * tempValue;
                saved = left[j - r] * tempValue;
            }
            derivatives[j] = saved;
        }
        return derivatives;
    }
    getParameterRange() {
        const uMin = this.knots[this.degree];
        const uMax = this.knots[this.knots.length - this.degree - 1];
        return [
            uMin,
            uMax
        ];
    }
    rationalBezierBasis(u, span, p, knotVector) {
        const basis = [];
        const left = [];
        const right = [];
        for (let i = 0; i <= p; i++) {
            left[i] = 1;
            right[i] = 1;
        }
        for (let j = 1; j <= p; j++) {
            for (let i = 0; i < j; i++) {
                left[j] = left[j] * (u - knotVector[span + 1 - j + i]) / (knotVector[span + i + 1] - knotVector[span + 1 - j + i]);
                right[j] = right[j] * (knotVector[span + j + i] - u) / (knotVector[span + j + 1] - knotVector[span + i + 1]);
            }
        }
        for (let j = 0; j <= p; j++) {
            basis[j] = left[j] * right[p - j];
        }
        return basis;
    }
    calculateDepthToPlane(point, plane) {
        const normal = plane.normal;
        const pointOnPlane = plane.point;
        const d = normal.x * pointOnPlane.x + normal.y * pointOnPlane.y + normal.z * pointOnPlane.z;
        return d - (normal.x * point.x + normal.y * point.y + normal.z * point.z);
    }
    getDenominator(u, span) {
        const n = this.controlPoints.length - 1;
        const p = this.degree;
        const basisDerivatives = this.basisFunctionDerivatives(span, u, p, this.knots);
        let denominator = 0;
        for (let i = 0; i <= n; i++) {
            denominator += basisDerivatives[0][i] * this.weights[i];
        }
        return denominator;
    }
    calculateProjection(point, plane) {
        const normal = plane.normal;
        const d = normal.x * plane.point.x + normal.y * plane.point.y + normal.z * plane.point.z;
        const t = (d - (normal.x * point.x + normal.y * point.y + normal.z * point.z)) / (normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        return {
            x: point.x + t * normal.x,
            y: point.y + t * normal.y,
            z: point.z + t * normal.z
        };
    }
    intersectWithCurve(curve) {
        const intersectionPoints = [];
        const samplingSteps = 100;
        const deltaT = 1 / samplingSteps;
        for (let t1 = 0; t1 <= 1; t1 += deltaT) {
            const point1 = this.evaluate(t1);
            for (let t2 = 0; t2 <= 1; t2 += deltaT) {
                const point2 = curve.evaluate(t2);
                const distance = point1.distanceTo(point2);
                if (distance <= defaultTolerance) {
                    const intersectionPoint = new Point(point1.x, point1.y, point1.z);
                    if (!intersectionPoints.some(pt => pt.equals(intersectionPoint))) {
                        intersectionPoints.push(intersectionPoint);
                    }
                }
            }
        }
        return intersectionPoints;
    }
}
class NURBSSurface {
    constructor(degreeU, degreeV, controlPoints, weights, knotsU, knotsV, isTrimmed) {
        this.degreeU = degreeU;
        this.degreeV = degreeV;
        this.controlPoints = controlPoints;
        this.weights = weights;
        this.knotsU = knotsU;
        this.knotsV = knotsV;
        this.isTrimmed = isTrimmed;
    }
    evaluate(u, v) {
        const n = this.controlPoints.length - 1;
        const m = this.controlPoints[0].length - 1;
        const p = this.degreeU;
        const q = this.degreeV;
        let result = new Point(0, 0, 0);
        let denominator = 0;
        for (let i = 0; i <= n; i++) {
            for (let j = 0; j <= m; j++) {
                const spanU = this.findSpan(n, p, u);
                const spanV = this.findSpan(m, q, v);
                const basisU = this.basisFunction(spanU, u, p, this.knotsU);
                const basisV = this.basisFunction(spanV, v, q, this.knotsV);
                const weight = this.weights[i][j];
                denominator += weight * basisU * basisV;
                const point = this.controlPoints[i][j];
                result.x += point.x * weight * basisU * basisV;
                result.y += point.y * weight * basisU * basisV;
                result.z += point.z * weight * basisU * basisV;
            }
        }
        if (denominator !== 0) {
            result.x /= denominator;
            result.y /= denominator;
            result.z /= denominator;
        }
        return result;
    }
    normal(u, v) {
        const spanU = this.findSpan(this.controlPoints.length - 1, this.degreeU, u);
        const spanV = this.findSpan(this.controlPoints[0].length - 1, this.degreeV, v);
        const basisU = this.basisFunction(spanU, u, this.degreeU, this.knotsU);
        const basisV = this.basisFunction(spanV, v, this.degreeV, this.knotsV);
        const tangentU = this.calculateTangent(spanU, u, basisU, this.degreeU, this.knotsU);
        const tangentV = this.calculateTangent(spanV, v, basisV, this.degreeV, this.knotsV);
        const normal = this.crossProduct(tangentU, tangentV);
        return normal.normalize();
    }
    splitU(u) {
        const n = this.controlPoints.length - 1;
        const newControlPoints = [];
        const newWeights = [];
        for (let i = 0; i <= n; i++) {
            const span = this.findSpan(n, this.degreeU, u);
            const basisU = this.basisFunction(span, u, this.degreeU, this.knotsU);
            const newRow = [];
            const newWeightsRow = [];
            for (let j = 0; j <= this.degreeV; j++) {
                const point = this.controlPoints[i][j];
                newRow.push(point);
                newWeightsRow.push(this.weights[i][j] * basisU);
            }
            newControlPoints.push(newRow);
            newWeights.push(newWeightsRow);
        }
        const newKnotsU = this.knotsU.slice();
        return new NURBSSurface(this.degreeU, this.degreeV, newControlPoints, newWeights, newKnotsU, this.knotsV, this.isTrimmed);
    }
    splitV(v) {
        const n = this.controlPoints.length - 1;
        const newControlPoints = [];
        const newWeights = [];
        for (let i = 0; i <= n; i++) {
            const span = this.findSpan(this.controlPoints[0].length - 1, this.degreeV, v);
            const basisV = this.basisFunction(span, v, this.degreeV, this.knotsV);
            const newRow = [];
            const newWeightsRow = [];
            for (let j = 0; j <= this.degreeU; j++) {
                const point = this.controlPoints[i][j];
                newRow.push(point);
                newWeightsRow.push(this.weights[i][j] * basisV);
            }
            newControlPoints.push(newRow);
            newWeights.push(newWeightsRow);
        }
        const newKnotsV = this.knotsV.slice();
        return new NURBSSurface(this.degreeU, this.degreeV, newControlPoints, newWeights, this.knotsU, newKnotsV, this.isTrimmed);
    }
    transform(matrix) {
        const transformedControlPoints = this.controlPoints.map(row => row.map(point => {
            const [x, y, z] = [
                point.x,
                point.y,
                point.z
            ];
            return new Point(matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3], matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3], matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3]);
        }));
        this.controlPoints = transformedControlPoints;
    }
    isPointOnSurface(point, tolerance = defaultTolerance) {
        const uMin = this.knotsU[this.degreeU];
        const uMax = this.knotsU[this.knotsU.length - this.degreeU - 1];
        const vMin = this.knotsV[this.degreeV];
        const vMax = this.knotsV[this.knotsV.length - this.degreeV - 1];
        if (point.x < uMin || point.x > uMax || point.y < vMin || point.y > vMax) {
            return false;
        }
        const surfacePoint = this.evaluate(point.x, point.y);
        return surfacePoint.distanceTo(point) <= tolerance;
    }
    refineKnots(direction) {
        if (direction === 'U') {
            const newKnotsU = this._refineKnotVector(this.knotsU, this.degreeU);
            this.knotsU = newKnotsU;
        } else if (direction === 'V') {
            const newKnotsV = this._refineKnotVector(this.knotsV, this.degreeV);
            this.knotsV = newKnotsV;
        } else {
            throw new Error('Direction must be either "U" or "V".');
        }
    }
    getBoundingBox() {
        const minX = Math.min(...this.controlPoints.flat().map(point => point.x));
        const minY = Math.min(...this.controlPoints.flat().map(point => point.y));
        const minZ = Math.min(...this.controlPoints.flat().map(point => point.z));
        const maxX = Math.max(...this.controlPoints.flat().map(point => point.x));
        const maxY = Math.max(...this.controlPoints.flat().map(point => point.y));
        const maxZ = Math.max(...this.controlPoints.flat().map(point => point.z));
        return {
            min: new Point(minX, minY, minZ),
            max: new Point(maxX, maxY, maxZ)
        };
    }
    getArea(precision) {
        let area = 0;
        const uMin = this.knotsU[this.degreeU];
        const uMax = this.knotsU[this.knotsU.length - this.degreeU - 1];
        const vMin = this.knotsV[this.degreeV];
        const vMax = this.knotsV[this.knotsV.length - this.degreeV - 1];
        const uSteps = Math.ceil((uMax - uMin) / precision);
        const vSteps = Math.ceil((vMax - vMin) / precision);
        for (let i = 0; i < uSteps; i++) {
            for (let j = 0; j < vSteps; j++) {
                const u1 = uMin + i / uSteps * (uMax - uMin);
                const v1 = vMin + j / vSteps * (vMax - vMin);
                const u2 = uMin + (i + 1) / uSteps * (uMax - uMin);
                const v2 = vMin + (j + 1) / vSteps * (vMax - vMin);
                const p1 = this.evaluate(u1, v1);
                const p2 = this.evaluate(u1, v2);
                const p3 = this.evaluate(u2, v2);
                const p4 = this.evaluate(u2, v1);
                area += this.calculateTriangleArea(p1, p2, p3);
                area += this.calculateTriangleArea(p1, p3, p4);
            }
        }
        return area;
    }
    projectToSurface(surface) {
        const projectedControlPoints = this.controlPoints.map(row => {
            return row.map(point => {
                // Calculate projection of each control point onto the target surface
                const targetPoint = surface.evaluate(point.x, point.y);
                // Assuming x and y are parametric coordinates of the surface
                return new Point(targetPoint.x, targetPoint.y, point.z);
            });
        });
        // Keeping the original z-coordinate
        this.controlPoints = projectedControlPoints;
    }
    getIsoCurve(direction, value) {
        const isoCurveControlPoints = [];
        const n = this.controlPoints.length - 1;
        const m = this.controlPoints[0].length - 1;
        if (direction === 'U') {
            const spanU = this.findSpan(n, this.degreeU, value);
            for (let j = 0; j <= m; j++) {
                const basisU = this.basisFunction(spanU, value, this.degreeU, this.knotsU);
                isoCurveControlPoints.push(new Point(this.controlPoints[spanU][j].x * basisU, this.controlPoints[spanU][j].y * basisU, this.controlPoints[spanU][j].z * basisU));
            }
        } else if (direction === 'V') {
            const spanV = this.findSpan(m, this.degreeV, value);
            for (let i = 0; i <= n; i++) {
                const basisV = this.basisFunction(spanV, value, this.degreeV, this.knotsV);
                isoCurveControlPoints.push(new Point(this.controlPoints[i][spanV].x * basisV, this.controlPoints[i][spanV].y * basisV, this.controlPoints[i][spanV].z * basisV));
            }
        } else {
            throw new Error('Direction must be either "U" or "V".');
        }
        return new NURBSCurve(this.degreeV, isoCurveControlPoints, Array(isoCurveControlPoints.length).fill(1), [], false);
    }
    trim(boundaries) {
        const trimmedControlPoints = this.controlPoints.map(row => row.map(point => point.clone()));
        const trimmedWeights = this.weights.map(row => row.slice());
        boundaries.forEach(boundary => {
            const boundaryPts = boundary.controlPoints;
            for (let i = 0; i < trimmedControlPoints.length; i++) {
                for (let j = 0; j < trimmedControlPoints[i].length; j++) {
                    const pt = trimmedControlPoints[i][j];
                    if (this.isPointInBoundary(pt, boundaryPts)) {
                        trimmedControlPoints[i][j] = new Point(0, 0, 0);
                        // Or any other clear representation of a trimmed point
                        trimmedWeights[i][j] = 0;
                    }
                }
            }
        });
        // Set the weight to zero
        this.controlPoints = trimmedControlPoints;
        this.weights = trimmedWeights;
    }
    calculateTriangleArea(p1, p2, p3) {
        const a = p2.distanceTo(p1);
        const b = p3.distanceTo(p2);
        const c = p1.distanceTo(p3);
        const s = (a + b + c) / 2;
        return Math.sqrt(s * (s - a) * (s - b) * (s - c));
    }
    intersectWireWithFace(wire) {
        const intersectionPoints = [];
        this.controlPoints.forEach(row => {
            row.forEach(point => {
                const wireStart = wire.startVertex.point;
                const wireEnd = wire.endVertex.point;
                const t = this.calculateIntersection(point, wireStart, wireEnd);
                if (t >= 0 && t <= 1) {
                    intersectionPoints.push(point);
                }
            });
        });
        return intersectionPoints;
    }
    computeArea(precision) {
        let area = 0;
        const uMin = this.knotsU[this.degreeU];
        const uMax = this.knotsU[this.knotsU.length - this.degreeU - 1];
        const vMin = this.knotsV[this.degreeV];
        const vMax = this.knotsV[this.knotsV.length - this.degreeV - 1];
        const uSteps = Math.ceil((uMax - uMin) / precision);
        const vSteps = Math.ceil((vMax - vMin) / precision);
        for (let i = 0; i < uSteps; i++) {
            for (let j = 0; j < vSteps; j++) {
                const u1 = uMin + i / uSteps * (uMax - uMin);
                const v1 = vMin + j / vSteps * (vMax - vMin);
                const u2 = uMin + (i + 1) / uSteps * (uMax - uMin);
                const v2 = vMin + (j + 1) / vSteps * (vMax - vMin);
                const p1 = this.evaluate(u1, v1);
                const p2 = this.evaluate(u1, v2);
                const p3 = this.evaluate(u2, v2);
                const p4 = this.evaluate(u2, v1);
                area += this.calculateTriangleArea(p1, p2, p3);
                area += this.calculateTriangleArea(p1, p3, p4);
            }
        }
        return area;
    }
    calculateTangent(span, u, basis, degree, knotVector) {
        const tangent = new Point(0, 0, 0);
        const weights = this.weights;
        for (let i = 0; i <= span; i++) {
            if (i < degree) {
                const derivative = this.basisFunctionDerivatives(span, u, degree - 1, knotVector);
                tangent.x += derivative[0] * weights[i] * this.controlPoints[i][0].x;
                tangent.y += derivative[0] * weights[i] * this.controlPoints[i][0].y;
                tangent.z += derivative[0] * weights[i] * this.controlPoints[i][0].z;
            }
        }
        return tangent;
    }
    crossProduct(v1, v2) {
        return new Point(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }
    _refineKnotVector(knots, degree) {
        const newKnots = [];
        const n = knots.length - degree - 1;
        for (let i = 0; i <= n; i++) {
            const knotSpan = knots[i + 1] - knots[i];
            newKnots.push(knots[i]);
            const additionalKnotsCount = Math.floor(Math.max(1, knotSpan * 10));
            for (let j = 1; j <= additionalKnotsCount; j++) {
                const newKnot = knots[i] + j * (knotSpan / (additionalKnotsCount + 1));
                newKnots.push(newKnot);
            }
        }
        newKnots.push(knots[n + 1]);
        return newKnots;
    }
    isPointInBoundary(point, boundaryPts) {
        // Check if a point is inside the boundary defined by boundaryPts
        // This is a placeholder; actual implementation would depend on boundary representation
        const bbox = this.getBoundingBox();
        if (point.x < bbox.min.x || point.x > bbox.max.x || point.y < bbox.min.y || point.y > bbox.max.y || point.z < bbox.min.z || point.z > bbox.max.z) {
            return false;
        }
        // Additional point-in-polygon or point-in-surface logic can be added here
        return true;
    }
    length() {
        const uMin = this.knotsU[this.degreeU];
        const uMax = this.knotsU[this.knotsU.length - this.degreeU - 1];
        const vMin = this.knotsV[this.degreeV];
        const vMax = this.knotsV[this.knotsV.length - this.degreeV - 1];
        let totalLength = 0;
        const uSteps = 100;
        const vSteps = 100;
        for (let i = 0; i <= uSteps; i++) {
            for (let j = 0; j <= vSteps; j++) {
                const u = uMin + (uMax - uMin) * (i / uSteps);
                const v = vMin + (vMax - vMin) * (j / vSteps);
                const point = this.evaluate(u, v);
                if (i > 0) {
                    const prevU = uMin + (uMax - uMin) * ((i - 1) / uSteps);
                    const prevPoint = this.evaluate(prevU, v);
                    totalLength += point.distanceTo(prevPoint);
                }
                if (j > 0) {
                    const prevV = vMin + (vMax - vMin) * ((j - 1) / vSteps);
                    const prevPoint = this.evaluate(u, prevV);
                    totalLength += point.distanceTo(prevPoint);
                }
            }
        }
        return totalLength;
    }
    calculateIntersection(surfacePoint, wireStart, wireEnd) {
        const epsilon = 0.000001;
        // Tolerance for intersection
        const dX = wireEnd.x - wireStart.x;
        const dY = wireEnd.y - wireStart.y;
        const parametricLine = t => ({
            x: wireStart.x + t * dX,
            y: wireStart.y + t * dY
        });
        let closestT = null;
        let closestDistance = Infinity;
        for (let t = 0; t <= 1; t += epsilon) {
            const linePoint = parametricLine(t);
            const distance = surfacePoint.distanceTo(new Point(linePoint.x, linePoint.y, surfacePoint.z));
            if (distance < closestDistance) {
                closestDistance = distance;
                closestT = t;
            }
        }
        return closestT;
    }
}
class Edge {
    constructor(curve, startVertex, endVertex, id) {
        this.curve = curve;
        this.startVertex = startVertex;
        this.endVertex = endVertex;
        this.id = id;
    }
    length() {
        const start = this.startVertex.point;
        const end = this.endVertex.point;
        return start.distanceTo(end);
    }
    midpoint() {
        const start = this.startVertex.point;
        const end = this.endVertex.point;
        return new Point((start.x + end.x) / 2, (start.y + end.y) / 2, (start.z + end.z) / 2);
    }
    transform(matrix) {
        this.startVertex.transform(matrix);
        this.endVertex.transform(matrix);
        this.curve.transform(matrix);
    }
    split(t) {
        const startPoint = this.startVertex.point;
        const endPoint = this.endVertex.point;
        const newPoint = startPoint.lerp(endPoint, t);
        const newVertex = new Vertex(newPoint, null);
        const newCurve = this.curve.split(t)[1];
        // Assuming split returns two curves
        const newEdge = new Edge(newCurve, newVertex, this.endVertex, `${ this.id }_split`);
        this.endVertex = newVertex;
        return [
            this,
            newEdge
        ];
    }
    reverse() {
        const temp = this.startVertex;
        this.startVertex = this.endVertex;
        this.endVertex = temp;
    }
    isClosed() {
        return this.startVertex.equals(this.endVertex);
    }
    getTangent(t) {
        return this.curve.tangent(t);
    }
    subdivide(n) {
        const segments = [];
        const step = 1 / n;
        for (let i = 0; i < n; i++) {
            const t = i * step;
            const newEdge = new Edge(this.curve.split(t)[0], this.startVertex, this.curve.split(t)[1].startVertex, `${ this.id }_${ i }`);
            segments.push(newEdge);
        }
        return segments;
    }
    validate() {
        const start = this.startVertex.point;
        const end = this.endVertex.point;
        return !start.equals(end);
    }
}
class Wire {
    constructor(edges, isClosed, id) {
        this.edges = edges;
        this.isClosed = isClosed;
        this.id = id;
    }
    validate() {
        if (this.edges.length < 2) {
            return false;
        }
        // A wire must have at least two edges to be valid
        const startVertex = this.edges[0].startVertex;
        const endVertex = this.edges[this.edges.length - 1].endVertex;
        // For closed wires, the start and end vertices must be the same
        if (this.isClosed && !startVertex.equals(endVertex)) {
            return false;
        }
        // Check the connectivity of the edges
        for (let i = 0; i < this.edges.length - 1; i++) {
            if (!this.edges[i].endVertex.equals(this.edges[i + 1].startVertex)) {
                return false;
            }
        }
        // Edges must connect properly
        // Check if all edges are valid
        return this.edges.every(edge => edge.validate());
    }
    computeBoundingBox() {
        const minX = Math.min(...this.edges.map(edge => Math.min(edge.startVertex.point.x, edge.endVertex.point.x)));
        const minY = Math.min(...this.edges.map(edge => Math.min(edge.startVertex.point.y, edge.endVertex.point.y)));
        const minZ = Math.min(...this.edges.map(edge => Math.min(edge.startVertex.point.z, edge.endVertex.point.z)));
        const maxX = Math.max(...this.edges.map(edge => Math.max(edge.startVertex.point.x, edge.endVertex.point.x)));
        const maxY = Math.max(...this.edges.map(edge => Math.max(edge.startVertex.point.y, edge.endVertex.point.y)));
        const maxZ = Math.max(...this.edges.map(edge => Math.max(edge.startVertex.point.z, edge.endVertex.point.z)));
        return {
            min: new Point(minX, minY, minZ),
            max: new Point(maxX, maxY, maxZ)
        };
    }
    transform(matrix) {
        this.edges.forEach(edge => {
            edge.transform(matrix);
        });
    }
    addEdge(edge) {
        this.edges.push(edge);
    }
    removeEdge(edge) {
        this.edges = this.edges.filter(e => e !== edge);
    }
    isPlanar(tolerance) {
        const edges = this.edges;
        if (edges.length < 3)
            return false;
        const normal = edges[0].getTangent(0).cross(edges[1].getTangent(0)).normalize();
        return edges.every(edge => {
            const edgeNormal = edge.getTangent(0).cross(edge.getTangent(1)).normalize();
            return normal.equals(edgeNormal, tolerance);
        });
    }
    reverse() {
        this.edges.reverse();
        if (this.isClosed) {
            const startVertex = this.edges[0].startVertex;
            const endVertex = this.edges[this.edges.length - 1].endVertex;
            this.edges[0].startVertex = endVertex;
            this.edges[this.edges.length - 1].endVertex = startVertex;
        }
    }
    close() {
        if (!this.isClosed) {
            const firstVertex = this.edges[0].startVertex;
            const lastEdge = this.edges[this.edges.length - 1];
            const lastVertex = lastEdge.endVertex;
            if (!firstVertex.equals(lastVertex)) {
                const closingCurve = new NURBSCurve(1, [
                    lastVertex.point,
                    firstVertex.point
                ], [
                    1,
                    1
                ], [], true);
                const closingEdge = new Edge(closingCurve, lastVertex, firstVertex, `${ this.id }_closing_edge`);
                this.addEdge(closingEdge);
                this.isClosed = true;
            }
        }
    }
    getLength() {
        return this.edges.reduce((total, edge) => total + edge.length(), 0);
    }
    factoryMethod(edges, isClosed = false) {
        return new Wire(edges, isClosed, `Wire_${ Date.now() }`);
    }
    isClosed() {
        return this.isClosed;
    }
    hollow(thickness) {
        const newEdges = [];
        this.edges.forEach(edge => {
            const start = edge.startVertex.point;
            const end = edge.endVertex.point;
            const direction = new Point(end.x - start.x, end.y - start.y, end.z - start.z);
            const length = Math.sqrt(direction.x ** 2 + direction.y ** 2 + direction.z ** 2);
            const offset = new Point(direction.x / length * thickness, direction.y / length * thickness, direction.z / length * thickness);
            const newStartVertex = new Vertex(new Point(start.x + offset.x, start.y + offset.y, start.z + offset.z), null);
            const newEndVertex = new Vertex(new Point(end.x + offset.x, end.y + offset.y, end.z + offset.z), null);
            const newCurve = edge.curve;
            // clone curve or modify if necessary
            const newEdge = new Edge(newCurve, newStartVertex, newEndVertex, `${ edge.id }_hollow`);
            newEdges.push(newEdge);
        });
        return new Wire(newEdges, this.isClosed, `${ this.id }_hollow`);
    }
    // ... other methods
    repairGaps(tolerance) {
        const repairedEdges = [];
        for (let i = 0; i < this.edges.length; i++) {
            const currentEdge = this.edges[i];
            const nextEdge = this.edges[(i + 1) % this.edges.length];
            const distance = currentEdge.endVertex.point.distanceTo(nextEdge.startVertex.point);
            if (distance > tolerance) {
                const midPoint = currentEdge.midpoint();
                const newVertex = new Vertex(midPoint, null);
                repairedEdges.push(new Edge(currentEdge.curve, currentEdge.startVertex, newVertex, `${ currentEdge.id }_repair`));
                repairedEdges.push(new Edge(currentEdge.curve, newVertex, nextEdge.startVertex, `${ nextEdge.id }_repair`));
            } else {
                repairedEdges.push(currentEdge);
            }
        }
        this.edges = repairedEdges;
    }
    split(t) {
        const newEdges = [];
        const newVertices = [];
        let totalLength = this.getLength();
        const splitEdgeCount = Math.floor(totalLength);
        for (let i = 0; i < splitEdgeCount; i++) {
            const edge = this.edges[i];
            const startVertex = edge.startVertex;
            const endVertex = edge.endVertex;
            if (totalLength > 0) {
                const edgeLength = edge.length();
                const splitT = t * edgeLength;
                if (splitT < edgeLength) {
                    const splitPoint = startVertex.point.lerp(endVertex.point, splitT / edgeLength);
                    const newVertex = new Vertex(splitPoint, null);
                    newVertices.push(newVertex);
                    const newEdge = new Edge(edge.curve.split(splitT / edgeLength)[0], startVertex, newVertex, `${ edge.id }_split_${ i }`);
                    newEdges.push(newEdge);
                    edge.startVertex = newVertex;
                }
            }
        }
        // Update the start vertex for the next edge
        // Add the remaining edges
        newEdges.push(...this.edges.slice(splitEdgeCount));
        return new Wire(newEdges, this.isClosed, `${ this.id }_split`);
    }
}
class Face {
    constructor(surface, outerWire, innerWires, id) {
        this.surface = surface;
        this.outerWire = outerWire;
        this.innerWires = innerWires;
        this.id = id;
    }
    area() {
        return this.surface.getArea(0.01);
    }
    tessellate(precision) {
        const vertices = [];
        const indices = [];
        const [uMin, uMax] = this.surface.getParameterRange();
        const [vMin, vMax] = this.surface.getParameterRange('V');
        const uSteps = Math.ceil((uMax - uMin) / precision);
        const vSteps = Math.ceil((vMax - vMin) / precision);
        for (let i = 0; i <= uSteps; i++) {
            for (let j = 0; j <= vSteps; j++) {
                const u = uMin + (uMax - uMin) * (i / uSteps);
                const v = vMin + (vMax - vMin) * (j / vSteps);
                const point = this.surface.evaluate(u, v);
                vertices.push(point);
                if (i < uSteps && j < vSteps) {
                    const currentIndex = i * (vSteps + 1) + j;
                    const nextIndex = currentIndex + vSteps + 1;
                    indices.push(currentIndex, nextIndex, currentIndex + 1);
                    indices.push(nextIndex, nextIndex + 1, currentIndex + 1);
                }
            }
        }
        return {
            vertices,
            indices
        };
    }
    transform(matrix) {
        this.surface.transform(matrix);
        this.outerWire.transform(matrix);
        this.innerWires.forEach(innerWire => innerWire.transform(matrix));
    }
    containsPoint(point) {
        const isInsideOuterWire = this.outerWire.isPlanar(defaultTolerance) && this.isPointInBoundary(point, this.outerWire.computeBoundingBox());
        const isInsideInnerWires = this.innerWires.every(innerWire => innerWire.isPlanar(defaultTolerance) && this.isPointInBoundary(point, innerWire.computeBoundingBox()));
        return isInsideOuterWire && !isInsideInnerWires;
    }
    addInnerWire(wire) {
        if (wire && wire instanceof Wire) {
            this.innerWires.push(wire);
        } else {
            throw new Error('Invalid wire. Ensure it is an instance of Wire.');
        }
    }
    removeInnerWire(wire) {
        const index = this.innerWires.indexOf(wire);
        if (index !== -1) {
            this.innerWires.splice(index, 1);
        } else {
            throw new Error('Inner wire not found in the face.');
        }
    }
    getNormal(u, v) {
        const dU = this.surface.evaluate(u + defaultTolerance, v);
        const dV = this.surface.evaluate(u, v + defaultTolerance);
        const tangentU = new Point(dU.x - this.surface.evaluate(u, v).x, dU.y - this.surface.evaluate(u, v).y, dU.z - this.surface.evaluate(u, v).z);
        const tangentV = new Point(dV.x - this.surface.evaluate(u, v).x, dV.y - this.surface.evaluate(u, v).y, dV.z - this.surface.evaluate(u, v).z);
        return tangentU.cross(tangentV).normalize();
    }
    validate() {
        if (!this.outerWire.validate()) {
            return false;
        }
        // Outer wire must be valid
        for (const innerWire of this.innerWires) {
            if (!innerWire.validate()) {
                return false;
            }
        }
        // All inner wires must be valid
        return true;
    }
    splitCurve(curve) {
        const intersectionPoints = this.outerWire.intersectCurve(curve);
        const newWires = [];
        if (intersectionPoints.length > 0) {
            const splitCurves = this.outerWire.splitAt(intersectionPoints);
            newWires.push(...splitCurves);
        }
        const newFace = new Face(this.surface, newWires, this.innerWires, this.id);
        return newFace;
    }
    mergeAdjacentFaces(face) {
        // Check if both faces share a common wire
        const sharedEdges = this.outerWire.edges.filter(edge1 => face.outerWire.edges.some(edge2 => edge1.startVertex.equals(edge2.startVertex) && edge1.endVertex.equals(edge2.endVertex)));
        if (sharedEdges.length > 0) {
            // Merge outer wires
            const newEdges = [
                ...this.outerWire.edges,
                ...face.outerWire.edges.filter(edge => !sharedEdges.includes(edge))
            ];
            this.outerWire.edges = newEdges;
            // Merge inner wires
            const newInnerWires = [
                ...this.innerWires,
                ...face.innerWires
            ];
            this.innerWires = newInnerWires;
            // Optionally, you may want to set face as undefined if it is merged completely
            return true;
        }
        // Indicates successful merge
        return false;
    }
}
class Shell {
    constructor(faces, id) {
        this.faces = faces;
        this.id = id;
    }
    validate() {
        if (!this.faces.length) {
            return false;
        }
        // A shell must have faces
        return this.faces.every(face => face.validate());
    }
    computeVolume() {
        let volume = 0;
        this.faces.forEach(face => {
            volume += face.area();
        });
        // Simplified for demonstration purposes
        return volume;
    }
    transform(matrix) {
        this.faces.forEach(face => {
            face.transform(matrix);
        });
    }
    addFace(face) {
        this.faces.push(face);
    }
    removeFace(face) {
        this.faces = this.faces.filter(f => f !== face);
    }
    isClosed() {
        return this.faces.every(face => face.outerWire.isClosed);
    }
    merge(shell) {
        if (!shell || !(shell instanceof Shell)) {
            throw new Error('Invalid shell to merge.');
        }
        const mergedFaces = [...this.faces];
        shell.faces.forEach(face => {
            const existingFace = mergedFaces.find(f => f.id === face.id);
            if (!existingFace) {
                mergedFaces.push(face);
            } else {
                existingFace.mergeAdjacentFaces(face);
            }
        });
        return new Shell(mergedFaces, `${ this.id }_merged_${ shell.id }`);
    }
    repairGaps(tolerance) {
    }
    computeSurfaceArea() {
        let totalArea = 0;
        this.faces.forEach(face => {
            totalArea += face.area();
        });
        return totalArea;
    }
}
class Solid {
    constructor(shell, id) {
        this.shell = shell;
        this.id = id;
    }
    volume() {
        return this.shell.computeVolume();
    }
    surfaceArea() {
        return this.shell.computeSurfaceArea();
    }
    pointInside(point) {
        const bbox = this.getBoundingBox();
        const isInsideBoundingBox = point.x >= bbox.min.x && point.x <= bbox.max.x && point.y >= bbox.min.y && point.y <= bbox.max.y && point.z >= bbox.min.z && point.z <= bbox.max.z;
        if (!isInsideBoundingBox) {
            return false;
        }
        let insideCount = 0;
        const vertices = this.shell.faces.flatMap(face => face.outerWire.edges.map(edge => [
            edge.startVertex.point,
            edge.endVertex.point
        ]));
        const n = vertices.length;
        for (let i = 0, j = n - 1; i < n; j = i++) {
            const vi = vertices[i][0];
            const vj = vertices[j][0];
            if (vi.y > point.y !== vj.y > point.y && point.x < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y) + vi.x) {
                insideCount++;
            }
        }
        return insideCount % 2 === 1;
    }
    split(plane) {
    }
    getBoundingBox() {
        const bbox = {
            min: new Point(Infinity, Infinity, Infinity),
            max: new Point(-Infinity, -Infinity, -Infinity)
        };
        this.shell.faces.forEach(face => {
            const faceBBox = face.outerWire.computeBoundingBox();
            bbox.min.x = Math.min(bbox.min.x, faceBBox.min.x);
            bbox.min.y = Math.min(bbox.min.y, faceBBox.min.y);
            bbox.min.z = Math.min(bbox.min.z, faceBBox.min.z);
            bbox.max.x = Math.max(bbox.max.x, faceBBox.max.x);
            bbox.max.y = Math.max(bbox.max.y, faceBBox.max.y);
            bbox.max.z = Math.max(bbox.max.z, faceBBox.max.z);
        });
        return bbox;
    }
    intersectWithSolid(solid) {
    }
    hollow(thickness) {
        const newShell = new Shell([]);
        this.shell.faces.forEach(face => {
            const outerWire = face.outerWire.hollow(thickness);
            newShell.addFace(new Face(face.surface, outerWire, face.innerWires, face.id));
        });
        return new Solid(newShell);
    }
    validate() {
    }
    createSolidsFromIntersection(intersectionEdges) {
    }
    transform(matrix) {
        this.shell.faces.forEach(face => {
            face.transform(matrix);
        });
    }
}
class BooleanOperation {
    union(solid1, solid2) {
        const newShell = new Shell([]);
        const mergedFaces = [
            ...solid1.shell.faces,
            ...solid2.shell.faces
        ];
        const uniqueFaces = new Map();
        mergedFaces.forEach(face => {
            if (!uniqueFaces.has(face.id)) {
                uniqueFaces.set(face.id, face.clone());
            } else {
                uniqueFaces.get(face.id).mergeAdjacentFaces(face);
            }
        });
        uniqueFaces.forEach(face => newShell.addFace(face));
        return new Solid(newShell);
    }
    intersection(solid1, solid2) {
        const newShell = new Shell([]);
        return new Solid(newShell);
    }
    difference(solid1, solid2) {
        const newShell = new Shell([]);
        return new Solid(newShell);
    }
    split(solid, plane) {
        const intersectionTool = new IntersectionTool();
        const intersectionEdges = intersectionTool.intersectWithSolid(solid, plane);
        const newSolids = this.createSolidsFromIntersection(intersectionEdges);
        return newSolids;
    }
    mergeSolids(solids) {
        const newShell = new Shell([]);
        const uniqueFaces = new Map();
        solids.forEach(solid => {
            solid.shell.faces.forEach(face => {
                if (!uniqueFaces.has(face.id)) {
                    uniqueFaces.set(face.id, face);
                    newShell.addFace(face.clone());
                } else {
                    const existingFace = uniqueFaces.get(face.id);
                    existingFace.mergeAdjacentFaces(face);
                }
            });
        });
        return new Solid(newShell);
    }
    subtractMultiple(solid, solidsToSubtract) {
        let resultSolid = solid;
        solidsToSubtract.forEach(solidToSubtract => {
            const intersectionTool = new IntersectionTool();
            const intersectionEdges = intersectionTool.intersect(resultSolid, solidToSubtract);
            const newShell = new Shell([]);
            intersectionEdges.forEach(edge => {
                const face = new Face();
                newShell.addFace(face);
            });
            resultSolid = new Solid(newShell);
        });
        return resultSolid;
    }
    resolveCoplanarFaces(shell) {
        const coplanarFaces = new Map();
        shell.faces.forEach(face => {
            const key = face.surface.getNormal(0, 0).toArray().join(',');
            if (!coplanarFaces.has(key)) {
                coplanarFaces.set(key, []);
            }
            coplanarFaces.get(key).push(face);
        });
        coplanarFaces.forEach(faces => {
            if (faces.length > 1) {
                const mergedFace = faces.reduce((acc, face) => {
                    acc.mergeAdjacentFaces(face);
                    return acc;
                }, faces[0].clone());
                shell.removeFace(faces[0]);
                shell.addFace(mergedFace);
            }
        });
    }
    validateBooleanResult(solid) {
        if (!solid.shell || !solid.shell.faces.length) {
            throw new Error('Solid must have a valid shell with faces.');
        }
        const normals = new Set();
        solid.shell.faces.forEach(face => {
            const normal = face.surface.getNormal(0, 0);
            normals.add(normal.toArray().join(','));
        });
        if (normals.size > 1) {
            throw new Error('Solid has non-consistent normals; it may have self-intersections.');
        }
        const selfIntersections = this.findSelfIntersections(solid);
        if (selfIntersections.length) {
            throw new Error('Solid has self-intersections that are not resolved.');
        }
        return true;
    }
    simplifyBooleanResult(solid) {
        const uniqueFaces = new Map();
        const simplifiedShell = new Shell([]);
        solid.shell.faces.forEach(face => {
            const key = face.surface.getNormal(0, 0).toArray().join(',');
            if (!uniqueFaces.has(key)) {
                uniqueFaces.set(key, face.clone());
            } else {
                const existingFace = uniqueFaces.get(key);
                existingFace.mergeAdjacentFaces(face);
            }
        });
        uniqueFaces.forEach(face => simplifiedShell.addFace(face));
        return new Solid(simplifiedShell);
    }
    findSelfIntersections(solid) {
        const intersections = [];
        const faces = solid.shell.faces;
        for (let i = 0; i < faces.length; i++) {
            for (let j = i + 1; j < faces.length; j++) {
                const intersection = this.intersectFace(faces[i], faces[j]);
                if (intersection) {
                    intersections.push(intersection);
                }
            }
        }
        return intersections;
    }
}
class IntersectionTool {
    intersect(solid1, solid2) {
        const intersectionEdges = [];
        solid1.shell.faces.forEach(face1 => {
            solid2.shell.faces.forEach(face2 => {
                const intersection = this.intersectFace(face1, face2);
                if (intersection) {
                    intersectionEdges.push(...intersection);
                }
            });
        });
        return intersectionEdges;
    }
    classifyPoint(point, solid) {
        const vertices = solid.shell.faces.flatMap(face => face.outerWire.edges.map(edge => [
            edge.startVertex.point,
            edge.endVertex.point
        ]));
        let inside = false;
        const n = vertices.length;
        for (let i = 0, j = n - 1; i < n; j = i++) {
            const vi = vertices[i][0];
            const vj = vertices[j][0];
            if (vi.y > point.y !== vj.y > point.y && point.x < (vj.x - vi.x) * (point.y - vi.y) / (vj.y - vi.y) + vi.x) {
                inside = !inside;
            }
        }
        return inside ? 'inside' : 'outside';
    }
    intersectEdgeAndFace(edge, face) {
        return face.intersectWireWithFace(edge);
    }
    intersectShells(shell1, shell2) {
        const intersectionEdges = [];
        const faces1 = shell1.faces;
        const faces2 = shell2.faces;
        faces1.forEach(face1 => {
            faces2.forEach(face2 => {
                const intersections = this.intersectFace(face1, face2);
                if (intersections.length > 0) {
                    intersectionEdges.push(...intersections);
                }
            });
        });
        return intersectionEdges;
    }
    findOverlappingRegions(shell1, shell2) {
        const overlappingRegions = [];
        shell1.faces.forEach(face1 => {
            shell2.faces.forEach(face2 => {
                const intersections = this.intersectFace(face1, face2);
                if (intersections.length > 0) {
                    overlappingRegions.push({
                        face1,
                        face2,
                        intersections
                    });
                }
            });
        });
        return overlappingRegions;
    }
    intersectWireWithFace(wire, face) {
        const intersectionPoints = [];
        const edges = wire.edges;
        edges.forEach(edge => {
            const intersection = this.intersectEdgeWithFace(edge, face);
            intersection.forEach(point => intersectionPoints.push(point));
        });
        return intersectionPoints;
    }
    intersectFace(face1, face2) {
        const intersections = [];
        const edges1 = face1.outerWire.edges;
        const edges2 = face2.outerWire.edges;
        edges1.forEach(edge1 => {
            edges2.forEach(edge2 => {
                const intersectionPoints = this.intersectEdgeWithEdge(edge1, edge2);
                if (intersectionPoints.length > 0) {
                    intersections.push(...intersectionPoints);
                }
            });
        });
        return intersections;
    }
    intersectCurves(curve1, curve2) {
        const intersectionPoints = [];
        const steps = 100;
        // Define a number of steps for the parameterization
        const tStep = 1 / steps;
        for (let i = 0; i <= steps; i++) {
            const t1 = i * tStep;
            const point1 = curve1.evaluate(t1);
            for (let j = 0; j <= steps; j++) {
                const t2 = j * tStep;
                const point2 = curve2.evaluate(t2);
                const distance = point1.distanceTo(point2);
                if (distance <= defaultTolerance) {
                    intersectionPoints.push(point1);
                }
            }
        }
        return intersectionPoints;
    }
    intersectEdgeWithEdge(edge1, edge2) {
        const intersectionPoints = [];
        const start1 = edge1.startVertex.point;
        const end1 = edge1.endVertex.point;
        const start2 = edge2.startVertex.point;
        const end2 = edge2.endVertex.point;
        const t = this.calculateIntersection(start1, end1, start2, end2);
        if (t >= 0 && t <= 1) {
            const intersectionPoint = start1.lerp(end1, t);
            intersectionPoints.push(intersectionPoint);
        }
        return intersectionPoints;
    }
    calculateIntersection(start1, end1, start2, end2) {
        const x1 = start1.x, y1 = start1.y;
        const x2 = end1.x, y2 = end1.y;
        const x3 = start2.x, y3 = start2.y;
        const x4 = end2.x, y4 = end2.y;
        const denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
        if (denom === 0)
            return -1;
        // Lines are parallel
        const t = ((x3 - x1) * (y4 - y3) - (y3 - y1) * (x4 - x3)) / denom;
        const u = -((x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)) / denom;
        return t >= 0 && t <= 1 ? t : -1;
    }
    intersectEdgeWithFace(edge, face) {
        const intersectionPoints = [];
        const edgePoints = [
            edge.startVertex.point,
            edge.endVertex.point
        ];
        const outerWireEdges = face.outerWire.edges;
        outerWireEdges.forEach(outerEdge => {
            const startEdgePoint = edgePoints[0];
            const endEdgePoint = edgePoints[1];
            const startOuterEdgePoint = outerEdge.startVertex.point;
            const endOuterEdgePoint = outerEdge.endVertex.point;
            const t = this.calculateIntersection(startEdgePoint, endEdgePoint, startOuterEdgePoint, endOuterEdgePoint);
            if (t >= 0 && t <= 1) {
                const intersectionPoint = startEdgePoint.lerp(endEdgePoint, t);
                intersectionPoints.push(intersectionPoint);
            }
        });
        return intersectionPoints;
    }
}
class ShapeHealingTool {
    healGaps(wire, tolerance) {
        const repairedEdges = [];
        for (let i = 0; i < wire.edges.length; i++) {
            const currentEdge = wire.edges[i];
            const nextEdge = wire.edges[(i + 1) % wire.edges.length];
            const distance = currentEdge.endVertex.point.distanceTo(nextEdge.startVertex.point);
            if (distance > tolerance) {
                const midPoint = currentEdge.midpoint();
                const newEdge = new Edge(currentEdge.curve, currentEdge.startVertex, new Vertex(midPoint, null), `${ currentEdge.id }_repair`);
                repairedEdges.push(newEdge);
            }
        }
        wire.edges.push(...repairedEdges);
    }
    healIntersections(shell) {
        const overlaps = new Set();
        shell.faces.forEach(face => {
            const edges = face.outerWire.edges;
            edges.forEach(edge => {
                const key = `${ edge.startVertex.point.id }-${ edge.endVertex.point.id }`;
                if (overlaps.has(key)) {
                    overlaps.delete(key);
                } else {
                    overlaps.add(key);
                }
            });
        });
        overlaps.forEach(key => {
            const [startId, endId] = key.split('-');
            const startVertex = shell.faces.flatMap(face => face.outerWire.edges).find(edge => edge.startVertex.point.id === startId);
            const endVertex = shell.faces.flatMap(face => face.outerWire.edges).find(edge => edge.endVertex.point.id === endId);
            if (startVertex && endVertex) {
                const midPoint = startVertex.startVertex.point.lerp(endVertex.endVertex.point, 0.5);
                const newVertex = new Vertex(midPoint, null);
                const repairedEdge = new Edge(startVertex.curve, startVertex.startVertex, newVertex, `${ startVertex.id }_healed`);
                shell.faces.forEach(face => {
                    face.outerWire.addEdge(repairedEdge);
                });
            }
        });
    }
    simplifyTopology(entity) {
        const uniqueVertices = new Map();
        const edgesToRemove = new Set();
        entity.shell.faces.forEach(face => {
            face.outerWire.edges.forEach(edge => {
                const key = `${ edge.startVertex.point.id }-${ edge.endVertex.point.id }`;
                if (!uniqueVertices.has(key)) {
                    uniqueVertices.set(key, edge);
                } else {
                    edgesToRemove.add(key);
                }
            });
        });
        edgesToRemove.forEach(key => {
            const edge = uniqueVertices.get(key);
            const face = entity.shell.faces.find(f => f.outerWire.edges.includes(edge));
            if (face) {
                face.outerWire.removeEdge(edge);
            }
        });
        entity.shell.faces.forEach(face => {
            face.innerWires.forEach(innerWire => {
                innerWire.edges.forEach(edge => {
                    const key = `${ edge.startVertex.point.id }-${ edge.endVertex.point.id }`;
                    if (uniqueVertices.has(key)) {
                        innerWire.removeEdge(edge);
                    }
                });
            });
        });
    }
    removeDuplicateVertices(shell, tolerance) {
        const uniqueVertices = new Map();
        const verticesToRemove = [];
        shell.faces.forEach(face => {
            face.outerWire.edges.forEach(edge => {
                const startVertexKey = `${ edge.startVertex.point.x },${ edge.startVertex.point.y },${ edge.startVertex.point.z }`;
                const endVertexKey = `${ edge.endVertex.point.x },${ edge.endVertex.point.y },${ edge.endVertex.point.z }`;
                if (!uniqueVertices.has(startVertexKey)) {
                    uniqueVertices.set(startVertexKey, edge.startVertex);
                } else {
                    verticesToRemove.push(edge.startVertex);
                }
                if (!uniqueVertices.has(endVertexKey)) {
                    uniqueVertices.set(endVertexKey, edge.endVertex);
                } else {
                    verticesToRemove.push(edge.endVertex);
                }
            });
        });
        verticesToRemove.forEach(vertex => {
            shell.faces.forEach(face => {
                face.outerWire.edges.forEach(edge => {
                    if (edge.startVertex.equals(vertex, tolerance)) {
                        edge.startVertex = uniqueVertices.get(`${ vertex.point.x },${ vertex.point.y },${ vertex.point.z }`);
                    }
                    if (edge.endVertex.equals(vertex, tolerance)) {
                        edge.endVertex = uniqueVertices.get(`${ vertex.point.x },${ vertex.point.y },${ vertex.point.z }`);
                    }
                });
            });
        });
    }
    repairDegenerateFaces(face) {
        const vertices = face.outerWire.edges.flatMap(edge => [
            edge.startVertex,
            edge.endVertex
        ]);
        const uniqueVertices = new Map();
        const degenerateEdges = [];
        vertices.forEach(vertex => {
            const key = vertex.point.toArray().join(',');
            if (!uniqueVertices.has(key)) {
                uniqueVertices.set(key, vertex);
            } else {
                degenerateEdges.push(vertex);
            }
        });
        if (degenerateEdges.length > 0) {
            const repairedEdges = [];
            const uniqueVertexList = Array.from(uniqueVertices.values());
            for (let edge of face.outerWire.edges) {
                const newStartVertex = uniqueVertices.get(edge.startVertex.point.toArray().join(','));
                const newEndVertex = uniqueVertices.get(edge.endVertex.point.toArray().join(','));
                if (!repairedEdges.some(e => e.startVertex.equals(newStartVertex) && e.endVertex.equals(newEndVertex))) {
                    repairedEdges.push(new Edge(edge.curve, newStartVertex, newEndVertex, edge.id));
                }
            }
            face.outerWire.edges = repairedEdges;
        }
    }
    optimizeShell(shell) {
        const uniqueVertices = new Map();
        const edgesToRemove = new Set();
        shell.faces.forEach(face => {
            face.outerWire.edges.forEach(edge => {
                const key = `${ edge.startVertex.point.x },${ edge.startVertex.point.y },${ edge.startVertex.point.z }-${ edge.endVertex.point.x },${ edge.endVertex.point.y },${ edge.endVertex.point.z }`;
                if (!uniqueVertices.has(key)) {
                    uniqueVertices.set(key, edge);
                } else {
                    edgesToRemove.add(key);
                }
            });
        });
        edgesToRemove.forEach(key => {
            const edge = uniqueVertices.get(key);
            const face = shell.faces.find(f => f.outerWire.edges.includes(edge));
            if (face) {
                face.outerWire.removeEdge(edge);
            }
        });
        shell.faces.forEach(face => {
            face.innerWires.forEach(innerWire => {
                innerWire.edges.forEach(edge => {
                    const key = `${ edge.startVertex.point.x },${ edge.startVertex.point.y },${ edge.startVertex.point.z }-${ edge.endVertex.point.x },${ edge.endVertex.point.y },${ edge.endVertex.point.z }`;
                    if (uniqueVertices.has(key)) {
                        innerWire.removeEdge(edge);
                    }
                });
            });
        });
    }
    fixNormals(shell) {
        shell.faces.forEach(face => {
            const normal = face.surface.getNormal(0, 0);
            face.outerWire.edges.forEach(edge => {
                const startVertexNormal = edge.startVertex.point.normal;
                const endVertexNormal = edge.endVertex.point.normal;
                if (!startVertexNormal.equals(normal)) {
                    edge.startVertex.point.normal = normal;
                }
                if (!endVertexNormal.equals(normal)) {
                    edge.endVertex.point.normal = normal;
                }
            });
        });
    }
    validateShapeIntegrity(entity) {
        if (!entity.shell || !entity.shell.faces.length) {
            throw new Error('Entity must have a valid shell with faces.');
        }
        const edgeMap = new Map();
        entity.shell.faces.forEach(face => {
            face.outerWire.edges.forEach(edge => {
                const key = `${ edge.startVertex.point.id }-${ edge.endVertex.point.id }`;
                if (!edgeMap.has(key)) {
                    edgeMap.set(key, edge);
                } else {
                    throw new Error('Entity contains duplicate edges.');
                }
            });
        });
        entity.shell.faces.forEach(face => {
            if (face.outerWire.isClosed && !face.outerWire.validate()) {
                throw new Error('Face outer wire is not valid.');
            }
            face.innerWires.forEach(innerWire => {
                if (!innerWire.validate()) {
                    throw new Error('Inner wire is not valid.');
                }
            });
        });
    }
}
class TessellationTool {
    tessellate(face, precision) {
        const vertices = [];
        const indices = [];
        const [uMin, uMax] = face.surface.getParameterRange('U');
        const [vMin, vMax] = face.surface.getParameterRange('V');
        const uSteps = Math.ceil((uMax - uMin) / precision);
        const vSteps = Math.ceil((vMax - vMin) / precision);
        for (let i = 0; i <= uSteps; i++) {
            for (let j = 0; j <= vSteps; j++) {
                const u = uMin + (uMax - uMin) * (i / uSteps);
                const v = vMin + (vMax - vMin) * (j / vSteps);
                const point = face.surface.evaluate(u, v);
                vertices.push(point);
                if (i < uSteps && j < vSteps) {
                    const currentIndex = i * (vSteps + 1) + j;
                    const nextIndex = currentIndex + vSteps + 1;
                    indices.push(currentIndex, nextIndex, currentIndex + 1);
                    indices.push(nextIndex, nextIndex + 1, currentIndex + 1);
                }
            }
        }
        return {
            vertices,
            indices
        };
    }
    optimizeMesh(mesh) {
        const uniqueVertices = new Map();
        const optimizedIndices = [];
        let indexCounter = 0;
        mesh.vertices.forEach(vertex => {
            const key = vertex.toArray().join(',');
            if (!uniqueVertices.has(key)) {
                uniqueVertices.set(key, indexCounter++);
            }
        });
        mesh.indices.forEach(originalIndex => {
            const vertex = mesh.vertices[originalIndex];
            const key = vertex.toArray().join(',');
            const optimizedIndex = uniqueVertices.get(key);
            optimizedIndices.push(optimizedIndex);
        });
        const optimizedVertices = Array.from(uniqueVertices.keys()).map(key => {
            const [x, y, z] = key.split(',').map(Number);
            return new Point(x, y, z);
        });
        return {
            vertices: optimizedVertices,
            indices: optimizedIndices
        };
    }
    subdivideMesh(mesh, maxEdgeLength) {
        const vertices = mesh.vertices.slice();
        const indices = mesh.indices.slice();
        const newIndices = [];
        const edgeMap = new Map();
        const addEdge = (v1, v2) => {
            const edgeKey = v1 < v2 ? `${ v1 }-${ v2 }` : `${ v2 }-${ v1 }`;
            if (!edgeMap.has(edgeKey)) {
                edgeMap.set(edgeKey, [
                    v1,
                    v2
                ]);
            }
        };
        for (let i = 0; i < indices.length; i += 3) {
            const [i0, i1, i2] = [
                indices[i],
                indices[i + 1],
                indices[i + 2]
            ];
            const p0 = vertices[i0];
            const p1 = vertices[i1];
            const p2 = vertices[i2];
            const edges = [
                {
                    v1: p0,
                    v2: p1
                },
                {
                    v1: p1,
                    v2: p2
                },
                {
                    v1: p2,
                    v2: p0
                }
            ];
            edges.forEach(edge => {
                const length = edge.v1.distanceTo(edge.v2);
                if (length > maxEdgeLength) {
                    const midPoint = edge.v1.lerp(edge.v2, 0.5);
                    const newIndex = vertices.length;
                    vertices.push(midPoint);
                    addEdge(indices[i], newIndex);
                    addEdge(newIndex, indices[i + 1]);
                    addEdge(newIndex, indices[i + 2]);
                } else {
                    addEdge(i0, i1);
                    addEdge(i1, i2);
                    addEdge(i2, i0);
                }
            });
        }
        edgeMap.forEach(edge => {
            newIndices.push(...edge);
        });
        return {
            vertices,
            indices: newIndices
        };
    }
    computeNormals(mesh) {
        const normals = new Array(mesh.vertices.length).fill(new Point(0, 0, 0));
        const faceCount = new Array(mesh.vertices.length).fill(0);
        for (let i = 0; i < mesh.indices.length; i += 3) {
            const i0 = mesh.indices[i];
            const i1 = mesh.indices[i + 1];
            const i2 = mesh.indices[i + 2];
            const v0 = mesh.vertices[i0];
            const v1 = mesh.vertices[i1];
            const v2 = mesh.vertices[i2];
            const edge1 = new Point(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
            const edge2 = new Point(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
            const normal = new Point(edge1.y * edge2.z - edge1.z * edge2.y, edge1.z * edge2.x - edge1.x * edge2.z, edge1.x * edge2.y - edge1.y * edge2.x);
            normals[i0].x += normal.x;
            normals[i0].y += normal.y;
            normals[i0].z += normal.z;
            faceCount[i0]++;
            normals[i1].x += normal.x;
            normals[i1].y += normal.y;
            normals[i1].z += normal.z;
            faceCount[i1]++;
            normals[i2].x += normal.x;
            normals[i2].y += normal.y;
            normals[i2].z += normal.z;
            faceCount[i2]++;
        }
        for (let i = 0; i < normals.length; i++) {
            if (faceCount[i] > 0) {
                normals[i].x /= faceCount[i];
                normals[i].y /= faceCount[i];
                normals[i].z /= faceCount[i];
                const length = Math.sqrt(normals[i].x ** 2 + normals[i].y ** 2 + normals[i].z ** 2);
                if (length > 0) {
                    normals[i].x /= length;
                    normals[i].y /= length;
                    normals[i].z /= length;
                }
            }
        }
        return normals;
    }
    generateUVMapping(face) {
        const uvMapping = [];
        const [uMin, uMax] = face.surface.getParameterRange('U');
        const [vMin, vMax] = face.surface.getParameterRange('V');
        const controlPoints = face.surface.controlPoints;
        const uSteps = controlPoints.length - 1;
        const vSteps = controlPoints[0].length - 1;
        for (let i = 0; i <= uSteps; i++) {
            for (let j = 0; j <= vSteps; j++) {
                const u = uMin + (uMax - uMin) * (i / uSteps);
                const v = vMin + (vMax - vMin) * (j / vSteps);
                uvMapping.push([
                    u,
                    v
                ]);
            }
        }
        return uvMapping;
    }
    validateMesh(mesh) {
        return true;
    }
}
class TransformationTool {
    translate(entity, vector) {
        const translationMatrix = [
            [
                1,
                0,
                0,
                vector.x
            ],
            [
                0,
                1,
                0,
                vector.y
            ],
            [
                0,
                0,
                1,
                vector.z
            ],
            [
                0,
                0,
                0,
                1
            ]
        ];
        entity.transform(translationMatrix);
    }
    rotate(entity, axis, angle) {
        const {x, y, z} = axis;
        const length = Math.sqrt(x * x + y * y + z * z);
        // Normalize the axis
        const nx = x / length;
        const ny = y / length;
        const nz = z / length;
        const cosTheta = Math.cos(angle);
        const sinTheta = Math.sin(angle);
        const oneMinusCos = 1 - cosTheta;
        const rotationMatrix = [
            [
                cosTheta + nx * nx * oneMinusCos,
                nx * ny * oneMinusCos - nz * sinTheta,
                nx * nz * oneMinusCos + ny * sinTheta,
                0
            ],
            [
                ny * nx * oneMinusCos + nz * sinTheta,
                cosTheta + ny * ny * oneMinusCos,
                ny * nz * oneMinusCos - nx * sinTheta,
                0
            ],
            [
                nz * nx * oneMinusCos - ny * sinTheta,
                nz * ny * oneMinusCos + nx * sinTheta,
                cosTheta + nz * nz * oneMinusCos,
                0
            ],
            [
                0,
                0,
                0,
                1
            ]
        ];
        entity.transform(rotationMatrix);
    }
    scale(entity, factor) {
        const scalingMatrix = [
            [
                factor,
                0,
                0,
                0
            ],
            [
                0,
                factor,
                0,
                0
            ],
            [
                0,
                0,
                factor,
                0
            ],
            [
                0,
                0,
                0,
                1
            ]
        ];
        entity.transform(scalingMatrix);
    }
    mirror(entity, plane) {
        const mirrorMatrix = [
            [
                1,
                0,
                0,
                0
            ],
            [
                0,
                1,
                0,
                0
            ],
            [
                0,
                0,
                1,
                0
            ],
            [
                0,
                0,
                0,
                1
            ]
        ];
        if (plane === 'xy') {
            mirrorMatrix[2][2] = -1;
        } else // Reflect across XY plane
        if (plane === 'yz') {
            mirrorMatrix[0][0] = -1;
        } else // Reflect across YZ plane
        if (plane === 'zx') {
            mirrorMatrix[1][1] = -1;
        } else
            // Reflect across ZX plane
            {
                throw new Error('Invalid plane. Use "xy", "yz", or "zx".');
            }
        entity.transform(mirrorMatrix);
    }
    combineTransformations(transformations) {
        let combinedMatrix = [
            [
                1,
                0,
                0,
                0
            ],
            [
                0,
                1,
                0,
                0
            ],
            [
                0,
                0,
                1,
                0
            ],
            [
                0,
                0,
                0,
                1
            ]
        ];
        transformations.forEach(transformation => {
            combinedMatrix = this.multiplyMatrices(combinedMatrix, transformation);
        });
        return combinedMatrix;
    }
    applyTransformationSequence(entity, sequence) {
        let combinedMatrix = [
            [
                1,
                0,
                0,
                0
            ],
            [
                0,
                1,
                0,
                0
            ],
            [
                0,
                0,
                1,
                0
            ],
            [
                0,
                0,
                0,
                1
            ]
        ];
        sequence.forEach(transformation => {
            combinedMatrix = this.multiplyMatrices(combinedMatrix, transformation);
        });
        entity.transform(combinedMatrix);
    }
    decomposeTransformation(matrix) {
        const translation = {
            x: matrix[0][3],
            y: matrix[1][3],
            z: matrix[2][3]
        };
        const scale = {
            x: Math.sqrt(matrix[0][0] ** 2 + matrix[0][1] ** 2 + matrix[0][2] ** 2),
            y: Math.sqrt(matrix[1][0] ** 2 + matrix[1][1] ** 2 + matrix[1][2] ** 2),
            z: Math.sqrt(matrix[2][0] ** 2 + matrix[2][1] ** 2 + matrix[2][2] ** 2)
        };
        const rotationMatrix = [
            [
                matrix[0][0] / scale.x,
                matrix[0][1] / scale.x,
                matrix[0][2] / scale.x
            ],
            [
                matrix[1][0] / scale.y,
                matrix[1][1] / scale.y,
                matrix[1][2] / scale.y
            ],
            [
                matrix[2][0] / scale.z,
                matrix[2][1] / scale.z,
                matrix[2][2] / scale.z
            ]
        ];
        const rotation = this.extractRotationFromMatrix(rotationMatrix);
        return {
            translation,
            rotation,
            scale
        };
    }
    removeDuplicateVertices(shell, tolerance) {
        const uniqueVertices = new Map();
        const verticesToRemove = [];
        shell.faces.forEach(face => {
            face.outerWire.edges.forEach(edge => {
                const startVertexKey = `${ edge.startVertex.point.x },${ edge.startVertex.point.y },${ edge.startVertex.point.z }`;
                const endVertexKey = `${ edge.endVertex.point.x },${ edge.endVertex.point.y },${ edge.endVertex.point.z }`;
                if (!uniqueVertices.has(startVertexKey)) {
                    uniqueVertices.set(startVertexKey, edge.startVertex);
                } else {
                    verticesToRemove.push(edge.startVertex);
                }
                if (!uniqueVertices.has(endVertexKey)) {
                    uniqueVertices.set(endVertexKey, edge.endVertex);
                } else {
                    verticesToRemove.push(edge.endVertex);
                }
            });
        });
        verticesToRemove.forEach(vertex => {
            shell.faces.forEach(face => {
                face.outerWire.edges.forEach(edge => {
                    if (edge.startVertex.equals(vertex, tolerance)) {
                        edge.startVertex = uniqueVertices.get(`${ vertex.point.x },${ vertex.point.y },${ vertex.point.z }`);
                    }
                    if (edge.endVertex.equals(vertex, tolerance)) {
                        edge.endVertex = uniqueVertices.get(`${ vertex.point.x },${ vertex.point.y },${ vertex.point.z }`);
                    }
                });
            });
        });
    }
    multiplyMatrices(a, b) {
        const result = [
            [
                0,
                0,
                0,
                0
            ],
            [
                0,
                0,
                0,
                0
            ],
            [
                0,
                0,
                0,
                0
            ],
            [
                0,
                0,
                0,
                0
            ]
        ];
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                result[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j] + a[i][3] * b[3][j];
            }
        }
        return result;
    }
    extractRotationFromMatrix(matrix) {
        const sy = Math.sqrt(matrix[0][0] ** 2 + matrix[1][0] ** 2);
        const singular = sy < 0.0001;
        // If sy is close to zero
        let x, y, z;
        if (!singular) {
            x = Math.atan2(matrix[2][1], matrix[2][2]);
            y = Math.atan2(-matrix[2][0], sy);
            z = Math.atan2(matrix[1][0], matrix[0][0]);
        } else {
            x = Math.atan2(-matrix[1][2], matrix[1][1]);
            y = Math.atan2(-matrix[2][0], sy);
            z = 0;
        }
        return {
            x,
            y,
            z
        };
    }
    fixNormals(shell) {
        shell.faces.forEach(face => {
            const normal = face.surface.getNormal(0, 0);
            face.outerWire.edges.forEach(edge => {
                const startVertexNormal = edge.startVertex.point.normal;
                const endVertexNormal = edge.endVertex.point.normal;
                if (!startVertexNormal.equals(normal)) {
                    edge.startVertex.point.normal = normal;
                }
                if (!endVertexNormal.equals(normal)) {
                    edge.endVertex.point.normal = normal;
                }
            });
        });
    }
}