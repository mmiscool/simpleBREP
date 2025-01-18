
export class Point3D {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    getX() {
        return this.x;
    }
    getY() {
        return this.y;
    }
    getZ() {
        return this.z;
    }
    setX(x) {
        this.x = x;
    }
    setY(y) {
        this.y = y;
    }
    setZ(z) {
        this.z = z;
    }
    toArray() {
        return [
            this.x,
            this.y,
            this.z
        ];
    }
    distanceTo(otherPoint) {
        const dx = this.x - otherPoint.x;
        const dy = this.y - otherPoint.y;
        const dz = this.z - otherPoint.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
    translate(vector) {
        return new Point3D(this.x + vector.getX(), this.y + vector.getY(), this.z + vector.getZ());
    }
    add(vector) {
        return new Point3D(this.x + vector.getX(), this.y + vector.getY(), this.z + vector.getZ());
    }
    sub(vector) {
        return new Point3D(this.x - vector.getX(), this.y - vector.getY(), this.z - vector.getZ());
    }
}
export class Vector3D {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    getX() {
        return this.x;
    }
    getY() {
        return this.y;
    }
    getZ() {
        return this.z;
    }
    setX(x) {
        this.x = x;
    }
    setY(y) {
        this.y = y;
    }
    setZ(z) {
        this.z = z;
    }
    toArray() {
        return [
            this.x,
            this.y,
            this.z
        ];
    }
    magnitude() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }
    normalize() {
        const mag = this.magnitude();
        if (mag === 0) {
            return new Vector3D(0, 0, 0);
        }
        return new Vector3D(this.x / mag, this.y / mag, this.z / mag);
    }
    dot(otherVector) {
        return this.x * otherVector.getX() + this.y * otherVector.getY() + this.z * otherVector.getZ();
    }
    cross(otherVector) {
        const x = this.y * otherVector.getZ() - this.z * otherVector.getY();
        const y = this.z * otherVector.getX() - this.x * otherVector.getZ();
        const z = this.x * otherVector.getY() - this.y * otherVector.getX();
        return new Vector3D(x, y, z);
    }
    add(otherVector) {
        return new Vector3D(this.x + otherVector.getX(), this.y + otherVector.getY(), this.z + otherVector.getZ());
    }
    sub(otherVector) {
        return new Vector3D(this.x - otherVector.getX(), this.y - otherVector.getY(), this.z - otherVector.getZ());
    }
    scale(scalar) {
        return new Vector3D(this.x * scalar, this.y * scalar, this.z * scalar);
    }
}
export class NurbsCurve {
    constructor(degree, controlPoints, knots, weights) {
        this.degree = degree;
        this.controlPoints = controlPoints;
        this.knots = knots;
        this.weights = weights;
    }
    getDegree() {
        return this.degree;
    }
    getControlPoints() {
        return this.controlPoints;
    }
    getKnots() {
        return this.knots;
    }
    getWeights() {
        return this.weights;
    }
    evaluate(u) {
        if (u < this.knots[0] || u > this.knots[this.knots.length - 1]) {
            return null;
        }
        return Utils.deBoor(this.degree, this.controlPoints, this.knots, this.weights, u);
    }
    evaluateDerivative(u, order) {
        if (u < this.knots[0] || u > this.knots[this.knots.length - 1]) {
            return null;
        }
        if (order < 0) {
            return null;
        }
        return Utils.deBoorDerivative(this.degree, this.controlPoints, this.knots, this.weights, u, order);
    }
    length() {
        let length = 0;
        const numSamples = 100;
        if (this.knots.length === 0) {
            return length;
        }
        const start = this.knots[0];
        const end = this.knots[this.knots.length - 1];
        const step = (end - start) / numSamples;
        let prevPoint = this.evaluate(start);
        for (let i = 1; i <= numSamples; i++) {
            const u = start + i * step;
            const currentPoint = this.evaluate(u);
            if (currentPoint && prevPoint) {
                length += currentPoint.distanceTo(prevPoint);
                prevPoint = currentPoint;
            }
        }
        return length;
    }
    tangeant(u) {
        const derivative = this.evaluateDerivative(u, 1);
        if (!derivative)
            return null;
        return derivative.normalize();
    }
    insertKnot(u, multiplicity) {
        if (multiplicity <= 0) {
            return;
        }
        const n = this.controlPoints.length - 1;
        const p = this.degree;
        const knots = this.knots;
        if (u < knots[0] || u > knots[knots.length - 1]) {
            return;
        }
        const span = Utils.findSpan(n, p, u, knots);
        let newKnots = [...knots];
        let newControlPoints = [...this.controlPoints];
        let newWeights = [...this.weights];
        for (let j = 0; j < multiplicity; j++) {
            newKnots.splice(span + 1, 0, u);
            const alpha = [];
            for (let i = 0; i <= p; i++) {
                const numerator = u - knots[span - p + i];
                const denominator = knots[span + 1 + i] - knots[span - p + i];
                if (denominator === 0) {
                    alpha[i] = 0;
                } else {
                    alpha[i] = numerator / denominator;
                }
            }
            let newControlPointsTemp = [];
            let newWeightsTemp = [];
            for (let i = 0; i <= n; i++) {
                if (i <= span - p) {
                    newControlPointsTemp.push(newControlPoints[i]);
                    newWeightsTemp.push(newWeights[i]);
                } else if (i > span) {
                    newControlPointsTemp.push(newControlPoints[i]);
                    newWeightsTemp.push(newWeights[i]);
                } else {
                    let cp = newControlPoints[i];
                    let w = newWeights[i];
                    const prevCp = newControlPoints[i - 1];
                    const prevW = newWeights[i - 1];
                    if (i === span - p + 1) {
                        const newPoint = new Point3D(alpha[1] * cp.x + (1 - alpha[1]) * prevCp.x, alpha[1] * cp.y + (1 - alpha[1]) * prevCp.y, alpha[1] * cp.z + (1 - alpha[1]) * prevCp.z);
                        const newWeight = alpha[1] * w + (1 - alpha[1]) * prevW;
                        newControlPointsTemp.push(newPoint);
                        newWeightsTemp.push(newWeight);
                        newControlPointsTemp.push(cp);
                        newWeightsTemp.push(w);
                    } else if (i > span - p + 1) {
                        const newPoint = new Point3D(alpha[i - (span - p)] * cp.x + (1 - alpha[i - (span - p)]) * newControlPointsTemp[newControlPointsTemp.length - 1].x, alpha[i - (span - p)] * cp.y + (1 - alpha[i - (span - p)]) * newControlPointsTemp[newControlPointsTemp.length - 1].y, alpha[i - (span - p)] * cp.z + (1 - alpha[i - (span - p)]) * newControlPointsTemp[newControlPointsTemp.length - 1].z);
                        const newWeight = alpha[i - (span - p)] * w + (1 - alpha[i - (span - p)]) * newWeightsTemp[newWeightsTemp.length - 1];
                        newControlPointsTemp.push(newPoint);
                        newWeightsTemp.push(newWeight);
                        if (i < span) {
                            newControlPointsTemp.push(cp);
                            newWeightsTemp.push(w);
                        }
                    } else {
                        newControlPointsTemp.push(cp);
                        newWeightsTemp.push(w);
                    }
                }
            }
            newControlPoints = newControlPointsTemp;
            newWeights = newWeightsTemp;
            this.controlPoints = newControlPoints;
            this.weights = newWeights;
        }
        this.knots = newKnots;
    }
    removeKnot(u, multiplicity) {
        if (multiplicity <= 0) {
            return;
        }
        const n = this.controlPoints.length - 1;
        const p = this.degree;
        const knots = this.knots;
        if (u < knots[0] || u > knots[knots.length - 1]) {
            return;
        }
        const span = Utils.findSpan(n, p, u, knots);
        let newKnots = [...knots];
        let newControlPoints = [...this.controlPoints];
        let newWeights = [...this.weights];
        let knotCount = 0;
        for (let i = 0; i < knots.length; i++) {
            if (knots[i] === u) {
                knotCount++;
            }
        }
        if (knotCount < multiplicity) {
            return;
        }
        for (let j = 0; j < multiplicity; j++) {
            newKnots.splice(span + 1, 1);
            let newControlPointsTemp = [];
            let newWeightsTemp = [];
            for (let i = 0; i <= n; i++) {
                if (i <= span - p) {
                    newControlPointsTemp.push(newControlPoints[i]);
                    newWeightsTemp.push(newWeights[i]);
                } else if (i > span + 1) {
                    newControlPointsTemp.push(newControlPoints[i]);
                    newWeightsTemp.push(newWeights[i]);
                } else {
                    const alpha = [];
                    for (let k = 0; k <= p; k++) {
                        const numerator = u - knots[span - p + k];
                        const denominator = knots[span + 1 + k] - knots[span - p + k];
                        if (denominator === 0) {
                            alpha[k] = 0;
                        } else {
                            alpha[k] = numerator / denominator;
                        }
                    }
                    if (span - p - 1 < 0) {
                        newControlPointsTemp.push(newControlPoints[i]);
                        newWeightsTemp.push(newWeights[i]);
                    } else {
                        const cp = newControlPoints[i];
                        const w = newWeights[i];
                        const nextCp = newControlPoints[i + 1];
                        const nextW = newWeights[i + 1];
                        const newPoint = new Point3D((cp.x - (1 - alpha[1]) * nextCp.x) / alpha[1], (cp.y - (1 - alpha[1]) * nextCp.y) / alpha[1], (cp.z - (1 - alpha[1]) * nextCp.z) / alpha[1]);
                        const newWeight = (w - (1 - alpha[1]) * nextW) / alpha[1];
                        newControlPointsTemp.push(newPoint);
                        newWeightsTemp.push(newWeight);
                    }
                }
            }
            newControlPoints = newControlPointsTemp;
            newWeights = newWeightsTemp;
            this.controlPoints = newControlPoints;
            this.weights = newWeights;
        }
        this.knots = newKnots;
    }
    splitCurve(u) {
        if (u < this.knots[0] || u > this.knots[this.knots.length - 1]) {
            return null;
        }
        const n = this.controlPoints.length - 1;
        const p = this.degree;
        const knots = this.knots;
        const span = Utils.findSpan(n, p, u, knots);
        let newKnots = [...knots];
        let newControlPoints = [...this.controlPoints];
        let newWeights = [...this.weights];
        let k = 0;
        for (let i = 0; i < knots.length; i++) {
            if (knots[i] === u) {
                k++;
            }
        }
        let multiplicity = p + 1 - k;
        if (multiplicity > 0) {
            this.insertKnot(u, multiplicity);
            newKnots = this.knots;
            newControlPoints = this.controlPoints;
            newWeights = this.weights;
        }
        const spanSplit = Utils.findSpan(newControlPoints.length - 1, p, u, newKnots);
        const curve1ControlPoints = newControlPoints.slice(0, spanSplit + 1);
        const curve1Weights = newWeights.slice(0, spanSplit + 1);
        const curve1Knots = newKnots.slice(0, spanSplit + p + 2);
        const curve2ControlPoints = newControlPoints.slice(spanSplit);
        const curve2Weights = newWeights.slice(spanSplit);
        const curve2Knots = newKnots.slice(spanSplit);
        const curve1 = new NurbsCurve(this.degree, curve1ControlPoints, curve1Knots, curve1Weights);
        const curve2 = new NurbsCurve(this.degree, curve2ControlPoints, curve2Knots, curve2Weights);
        return [
            curve1,
            curve2
        ];
    }
    reverse() {
        const reversedControlPoints = [...this.controlPoints].reverse();
        const reversedWeights = [...this.weights].reverse();
        const reversedKnots = [...this.knots].reverse().map(k => this.knots[this.knots.length - 1] - k);
        return new NurbsCurve(this.degree, reversedControlPoints, reversedKnots, reversedWeights);
    }
}
export class NurbsSurface {
    constructor(degreeU, degreeV, controlPoints, knotsU, knotsV, weights) {
        this.degreeU = degreeU;
        this.degreeV = degreeV;
        this.controlPoints = controlPoints;
        this.knotsU = knotsU;
        this.knotsV = knotsV;
        this.weights = weights;
    }
    getDegreeU() {
        return this.degreeU;
    }
    getDegreeV() {
        return this.degreeV;
    }
    getControlPoints() {
        return this.controlPoints;
    }
    getKnotsU() {
        return this.knotsU;
    }
    getKnotsV() {
        return this.knotsV;
    }
    getWeights() {
        return this.weights;
    }
    evaluate(u, v) {
        if (u < this.knotsU[0] || u > this.knotsU[this.knotsU.length - 1] || v < this.knotsV[0] || v > this.knotsV[this.knotsV.length - 1]) {
            return null;
        }
        const n = this.controlPoints.length - 1;
        const m = this.controlPoints[0].length - 1;
        const p = this.degreeU;
        const q = this.degreeV;
        const U = this.knotsU;
        const V = this.knotsV;
        const P = this.controlPoints;
        const W = this.weights;
        let Su = [];
        for (let i = 0; i <= m; i++) {
            let controlPoints = [];
            let weights = [];
            for (let j = 0; j <= n; j++) {
                controlPoints.push(P[j][i]);
                weights.push(W[j][i]);
            }
            Su.push(Utils.deBoor(p, controlPoints, U, weights, u));
        }
        let controlPoints = [...Su];
        let weights = Array(Su.length).fill(1);
        return Utils.deBoor(q, controlPoints, V, weights, v);
    }
    evaluatePartialDerivative(u, v, orderU, orderV) {
        if (u < this.knotsU[0] || u > this.knotsU[this.knotsU.length - 1] || v < this.knotsV[0] || v > this.knotsV[this.knotsV.length - 1]) {
            return null;
        }
        if (orderU < 0 || orderV < 0) {
            return null;
        }
        const n = this.controlPoints.length - 1;
        const m = this.controlPoints[0].length - 1;
        const p = this.degreeU;
        const q = this.degreeV;
        const U = this.knotsU;
        const V = this.knotsV;
        const P = this.controlPoints;
        const W = this.weights;
        let Su = [];
        for (let i = 0; i <= m; i++) {
            let controlPoints = [];
            let weights = [];
            for (let j = 0; j <= n; j++) {
                controlPoints.push(P[j][i]);
                weights.push(W[j][i]);
            }
            Su.push(Utils.deBoorDerivative(p, controlPoints, U, weights, u, orderU));
        }
        let controlPoints = [...Su];
        let weights = Array(Su.length).fill(1);
        const derivative = Utils.deBoorDerivative(q, controlPoints, V, weights, v, orderV);
        if (!derivative || isNaN(derivative.getX()) || isNaN(derivative.getY()) || isNaN(derivative.getZ())) {
            return null;
        }
        return derivative;
    }
    normal(u, v) {
        const derivativeU = this.evaluatePartialDerivative(u, v, 1, 0);
        const derivativeV = this.evaluatePartialDerivative(u, v, 0, 1);
        if (!derivativeU || !derivativeV) {
            return null;
        }
        const normalVector = derivativeU.cross(derivativeV);
        if (normalVector.magnitude() === 0) {
            return null;
        }
        if (isNaN(normalVector.getX()) || isNaN(normalVector.getY()) || isNaN(normalVector.getZ())) {
            return null;
        }
        return normalVector.normalize();
    }
    insertKnotU(u, multiplicity) {
        if (multiplicity <= 0) {
            return;
        }
        const n = this.controlPoints.length - 1;
        const p = this.degreeU;
        const knots = this.knotsU;
        if (u < knots[0] || u > knots[knots.length - 1]) {
            return;
        }
        const span = Utils.findSpan(n, p, u, knots);
        let newKnots = [...knots];
        let newControlPoints = this.controlPoints.map(row => [...row]);
        let newWeights = this.weights.map(row => [...row]);
        for (let j = 0; j < multiplicity; j++) {
            newKnots.splice(span + 1, 0, u);
            const alpha = [];
            for (let i = 0; i <= p; i++) {
                const numerator = u - knots[span - p + i];
                const denominator = knots[span + 1 + i] - knots[span - p + i];
                if (denominator === 0) {
                    alpha[i] = 0;
                } else {
                    alpha[i] = numerator / denominator;
                }
            }
            let newControlPointsTemp = [];
            let newWeightsTemp = [];
            for (let i = 0; i < newControlPoints.length; i++) {
                let rowControlPoints = [];
                let rowWeights = [];
                for (let k = 0; k < newControlPoints[i].length; k++) {
                    if (i <= span - p) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else if (i > span) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else {
                        const cp = newControlPoints[i][k];
                        const w = newWeights[i][k];
                        const prevCp = newControlPoints[i - 1][k];
                        const prevW = newWeights[i - 1][k];
                        const newPoint = new Point3D(alpha[i - (span - p)] * cp.x + (1 - alpha[i - (span - p)]) * prevCp.x, alpha[i - (span - p)] * cp.y + (1 - alpha[i - (span - p)]) * prevCp.y, alpha[i - (span - p)] * cp.z + (1 - alpha[i - (span - p)]) * prevCp.z);
                        const newWeight = alpha[i - (span - p)] * w + (1 - alpha[i - (span - p)]) * prevW;
                        rowControlPoints.push(newPoint);
                        rowWeights.push(newWeight);
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    }
                }
                newControlPointsTemp.push(rowControlPoints);
                newWeightsTemp.push(rowWeights);
            }
            newControlPoints = newControlPointsTemp;
            newWeights = newWeightsTemp;
            this.controlPoints = newControlPoints;
            this.weights = newWeights;
        }
        this.knotsU = newKnots;
    }
    insertKnotV(v, multiplicity) {
        if (multiplicity <= 0) {
            return;
        }
        const m = this.controlPoints[0].length - 1;
        const q = this.degreeV;
        const knots = this.knotsV;
        if (v < knots[0] || v > knots[knots.length - 1]) {
            return;
        }
        const span = Utils.findSpan(m, q, v, knots);
        let newKnots = [...knots];
        let newControlPoints = this.controlPoints.map(row => [...row]);
        let newWeights = this.weights.map(row => [...row]);
        for (let j = 0; j < multiplicity; j++) {
            newKnots.splice(span + 1, 0, v);
            const alpha = [];
            for (let i = 0; i <= q; i++) {
                const numerator = v - knots[span - q + i];
                const denominator = knots[span + 1 + i] - knots[span - q + i];
                if (denominator === 0) {
                    alpha[i] = 0;
                } else {
                    alpha[i] = numerator / denominator;
                }
            }
            let newControlPointsTemp = [];
            let newWeightsTemp = [];
            for (let i = 0; i < newControlPoints.length; i++) {
                let rowControlPoints = [];
                let rowWeights = [];
                for (let k = 0; k < newControlPoints[i].length; k++) {
                    if (k <= span - q) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else if (k > span) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else {
                        const cp = newControlPoints[i][k];
                        const w = newWeights[i][k];
                        const prevCp = newControlPoints[i][k - 1];
                        const prevW = newWeights[i][k - 1];
                        const newPoint = new Point3D(alpha[k - (span - q)] * cp.x + (1 - alpha[k - (span - q)]) * prevCp.x, alpha[k - (span - q)] * cp.y + (1 - alpha[k - (span - q)]) * prevCp.y, alpha[k - (span - q)] * cp.z + (1 - alpha[k - (span - q)]) * prevCp.z);
                        const newWeight = alpha[k - (span - q)] * w + (1 - alpha[k - (span - q)]) * prevW;
                        rowControlPoints.push(newPoint);
                        rowWeights.push(newWeight);
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    }
                }
                newControlPointsTemp.push(rowControlPoints);
                newWeightsTemp.push(rowWeights);
            }
            newControlPoints = newControlPointsTemp;
            newWeights = newWeightsTemp;
            this.controlPoints = newControlPoints;
            this.weights = newWeights;
        }
        this.knotsV = newKnots;
    }
    removeKnotU(u, multiplicity) {
        if (multiplicity <= 0) {
            return;
        }
        const n = this.controlPoints.length - 1;
        const p = this.degreeU;
        const knots = this.knotsU;
        if (u <= knots[0] || u >= knots[knots.length - 1]) {
            return;
        }
        const span = Utils.findSpan(n, p, u, knots);
        let newKnots = [...knots];
        let newControlPoints = this.controlPoints.map(row => [...row]);
        let newWeights = this.weights.map(row => [...row]);
        let knotCount = 0;
        for (let i = 0; i < knots.length; i++) {
            if (knots[i] === u) {
                knotCount++;
            }
        }
        if (knotCount < multiplicity) {
            return;
        }
        for (let j = 0; j < multiplicity; j++) {
            newKnots.splice(span + 1, 1);
            let newControlPointsTemp = [];
            let newWeightsTemp = [];
            for (let i = 0; i < newControlPoints.length; i++) {
                let rowControlPoints = [];
                let rowWeights = [];
                for (let k = 0; k < newControlPoints[i].length; k++) {
                    if (i <= span - p) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else if (i > span + 1) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else {
                        const alpha = [];
                        for (let l = 0; l <= p; l++) {
                            const numerator = u - knots[span - p + l];
                            const denominator = knots[span + 1 + l] - knots[span - p + l];
                            if (denominator === 0) {
                                alpha[l] = 0;
                            } else {
                                alpha[l] = numerator / denominator;
                            }
                        }
                        if (span - p - 1 < 0) {
                            rowControlPoints.push(newControlPoints[i][k]);
                            rowWeights.push(newWeights[i][k]);
                        } else {
                            const cp = newControlPoints[i][k];
                            const w = newWeights[i][k];
                            const nextCp = newControlPoints[i + 1][k];
                            const nextW = newWeights[i + 1][k];
                            const newPoint = new Point3D((cp.x - (1 - alpha[1]) * nextCp.x) / alpha[1], (cp.y - (1 - alpha[1]) * nextCp.y) / alpha[1], (cp.z - (1 - alpha[1]) * nextCp.z) / alpha[1]);
                            const newWeight = (w - (1 - alpha[1]) * nextW) / alpha[1];
                            rowControlPoints.push(newPoint);
                            rowWeights.push(newWeight);
                        }
                    }
                }
                newControlPointsTemp.push(rowControlPoints);
                newWeightsTemp.push(rowWeights);
            }
            newControlPoints = newControlPointsTemp;
            newWeights = newWeightsTemp;
            this.controlPoints = newControlPoints;
            this.weights = newWeights;
        }
        this.knotsU = newKnots;
    }
    removeKnotV(v, multiplicity) {
        if (multiplicity <= 0) {
            return;
        }
        const m = this.controlPoints[0].length - 1;
        const q = this.degreeV;
        const knots = this.knotsV;
        if (v <= knots[0] || v >= knots[knots.length - 1]) {
            return;
        }
        const span = Utils.findSpan(m, q, v, knots);
        let newKnots = [...knots];
        let newControlPoints = this.controlPoints.map(row => [...row]);
        let newWeights = this.weights.map(row => [...row]);
        let knotCount = 0;
        for (let i = 0; i < knots.length; i++) {
            if (knots[i] === v) {
                knotCount++;
            }
        }
        if (knotCount < multiplicity) {
            return;
        }
        for (let j = 0; j < multiplicity; j++) {
            newKnots.splice(span + 1, 1);
            let newControlPointsTemp = [];
            let newWeightsTemp = [];
            for (let i = 0; i < newControlPoints.length; i++) {
                let rowControlPoints = [];
                let rowWeights = [];
                for (let k = 0; k < newControlPoints[i].length; k++) {
                    if (k <= span - q) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else if (k > span + 1) {
                        rowControlPoints.push(newControlPoints[i][k]);
                        rowWeights.push(newWeights[i][k]);
                    } else {
                        const alpha = [];
                        for (let l = 0; l <= q; l++) {
                            const numerator = v - knots[span - q + l];
                            const denominator = knots[span + 1 + l] - knots[span - q + l];
                            if (denominator === 0) {
                                alpha[l] = 0;
                            } else {
                                alpha[l] = numerator / denominator;
                            }
                        }
                        if (span - q - 1 < 0) {
                            rowControlPoints.push(newControlPoints[i][k]);
                            rowWeights.push(newWeights[i][k]);
                        } else {
                            const cp = newControlPoints[i][k];
                            const w = newWeights[i][k];
                            const nextCp = newControlPoints[i][k + 1];
                            const nextW = newWeights[i][k + 1];
                            const newPoint = new Point3D((cp.x - (1 - alpha[1]) * nextCp.x) / alpha[1], (cp.y - (1 - alpha[1]) * nextCp.y) / alpha[1], (cp.z - (1 - alpha[1]) * nextCp.z) / alpha[1]);
                            const newWeight = (w - (1 - alpha[1]) * nextW) / alpha[1];
                            rowControlPoints.push(newPoint);
                            rowWeights.push(newWeight);
                        }
                    }
                }
                newControlPointsTemp.push(rowControlPoints);
                newWeightsTemp.push(rowWeights);
            }
            newControlPoints = newControlPointsTemp;
            newWeights = newWeightsTemp;
            this.controlPoints = newControlPoints;
            this.weights = newWeights;
        }
        this.knotsV = newKnots;
    }
    splitSurfaceU(u) {
        if (u < this.knotsU[0] || u > this.knotsU[this.knotsU.length - 1]) {
            return null;
        }
        const n = this.controlPoints.length - 1;
        const p = this.degreeU;
        const knots = this.knotsU;
        const m = this.controlPoints[0].length - 1;
        const q = this.degreeV;
        const knotsV = this.knotsV;
        const span = Utils.findSpan(n, p, u, knots);
        let newKnots = [...knots];
        let newControlPoints = this.controlPoints.map(row => [...row]);
        let newWeights = this.weights.map(row => [...row]);
        let k = 0;
        for (let i = 0; i < knots.length; i++) {
            if (knots[i] === u) {
                k++;
            }
        }
        let multiplicity = p + 1 - k;
        if (multiplicity > 0) {
            this.insertKnotU(u, multiplicity);
            newKnots = this.knotsU;
            newControlPoints = this.controlPoints;
            newWeights = this.weights;
        }
        const spanSplit = Utils.findSpan(newControlPoints.length - 1, p, u, newKnots);
        const surface1ControlPoints = newControlPoints.slice(0, spanSplit + 1);
        const surface1Weights = newWeights.slice(0, spanSplit + 1);
        const surface1Knots = newKnots.slice(0, spanSplit + p + 2);
        const surface2ControlPoints = newControlPoints.slice(spanSplit);
        const surface2Weights = newWeights.slice(spanSplit);
        const surface2Knots = newKnots.slice(spanSplit);
        const surface1 = new NurbsSurface(this.degreeU, this.degreeV, surface1ControlPoints, surface1Knots, knotsV, surface1Weights);
        const surface2 = new NurbsSurface(this.degreeU, this.degreeV, surface2ControlPoints, surface2Knots, knotsV, surface2Weights);
        return [
            surface1,
            surface2
        ];
    }
    splitSurfaceV(v) {
        if (v < this.knotsV[0] || v > this.knotsV[this.knotsV.length - 1]) {
            return null;
        }
        const n = this.controlPoints.length - 1;
        const p = this.degreeU;
        const knotsU = this.knotsU;
        const m = this.controlPoints[0].length - 1;
        const q = this.degreeV;
        const knots = this.knotsV;
        const span = Utils.findSpan(m, q, v, knots);
        let newKnots = [...knots];
        let newControlPoints = this.controlPoints.map(row => [...row]);
        let newWeights = this.weights.map(row => [...row]);
        let k = 0;
        for (let i = 0; i < knots.length; i++) {
            if (knots[i] === v) {
                k++;
            }
        }
        let multiplicity = q + 1 - k;
        if (multiplicity > 0) {
            this.insertKnotV(v, multiplicity);
            newKnots = this.knotsV;
            newControlPoints = this.controlPoints;
            newWeights = this.weights;
        }
        const spanSplit = Utils.findSpan(newControlPoints[0].length - 1, q, v, newKnots);
        const surface1ControlPoints = [];
        const surface1Weights = [];
        const surface2ControlPoints = [];
        const surface2Weights = [];
        for (let i = 0; i < newControlPoints.length; i++) {
            surface1ControlPoints.push(newControlPoints[i].slice(0, spanSplit + 1));
            surface1Weights.push(newWeights[i].slice(0, spanSplit + 1));
            surface2ControlPoints.push(newControlPoints[i].slice(spanSplit));
            surface2Weights.push(newWeights[i].slice(spanSplit));
        }
        const surface1Knots = newKnots.slice(0, spanSplit + q + 2);
        const surface2Knots = newKnots.slice(spanSplit);
        const surface1 = new NurbsSurface(this.degreeU, this.degreeV, surface1ControlPoints, knotsU, surface1Knots, surface1Weights);
        const surface2 = new NurbsSurface(this.degreeU, this.degreeV, surface2ControlPoints, knotsU, surface2Knots, surface2Weights);
        return [
            surface1,
            surface2
        ];
    }
}
export class Utils {
    basisFunction(i, degree, u, knots) {
        if (degree === 0) {
            if (u >= knots[i] && u < knots[i + 1]) {
                return 1;
            } else {
                return 0;
            }
        }
        if (knots[i + degree + 1] === knots[i]) {
            return 0;
        }
        if (knots[i + degree] === knots[i]) {
            return 0;
        }
        const firstPart = (u - knots[i]) / (knots[i + degree] - knots[i]) * this.basisFunction(i, degree - 1, u, knots);
        const secondPart = (knots[i + degree + 1] - u) / (knots[i + degree + 1] - knots[i + 1]) * this.basisFunction(i + 1, degree - 1, u, knots);
        return firstPart + secondPart;
    }
    basisFunctionDerivative(i, degree, u, knots, order) {
        if (order < 0) {
            return 0;
        }
        if (order === 0) {
            return this.basisFunction(i, degree, u, knots);
        }
        if (degree === 0) {
            return 0;
        }
        let firstPart = 0;
        let secondPart = 0;
        if (knots[i + degree] !== knots[i]) {
            firstPart = degree / (knots[i + degree] - knots[i]) * this.basisFunctionDerivative(i, degree - 1, u, knots, order - 1);
        }
        if (knots[i + degree + 1] !== knots[i + 1]) {
            secondPart = -degree / (knots[i + degree + 1] - knots[i + 1]) * this.basisFunctionDerivative(i + 1, degree - 1, u, knots, order - 1);
        }
        return firstPart + secondPart;
    }
    deBoor(degree, controlPoints, knots, weights, u) {
        const n = controlPoints.length - 1;
        const span = this.findSpan(n, degree, u, knots);
        const basisValues = [];
        for (let i = 0; i <= degree; i++) {
            basisValues.push(this.basisFunction(span - degree + i, degree, u, knots));
        }
        let point = new Point3D(0, 0, 0);
        let weightSum = 0;
        for (let i = 0; i <= degree; i++) {
            const weightedPoint = controlPoints[span - degree + i].scale(weights[span - degree + i] * basisValues[i]);
            point = point.add(weightedPoint);
            weightSum += weights[span - degree + i] * basisValues[i];
        }
        if (weightSum === 0) {
            return new Point3D(0, 0, 0);
        }
        return point.scale(1 / weightSum);
    }
    deBoorDerivative(degree, controlPoints, knots, weights, u, order) {
        const n = controlPoints.length - 1;
        const span = this.findSpan(n, degree, u, knots);
        let ders = [];
        for (let k = 0; k <= order; k++) {
            ders.push([]);
        }
        for (let i = 0; i <= degree; i++) {
            ders[0].push(controlPoints[span - degree + i]);
        }
        for (let k = 1; k <= order; k++) {
            for (let i = 0; i <= degree - k; i++) {
                let alpha = (u - knots[span - degree + i + k]) / (knots[span + 1 + i] - knots[span - degree + i + k]);
                if (isNaN(alpha))
                    alpha = 0;
                let point = ders[k - 1][i + 1].sub(ders[k - 1][i]).scale(degree - k + 1);
                ders[k].push(point.scale(1 / (knots[span + 1 + i] - knots[span - degree + i + k])));
            }
        }
        let vector = new Vector3D(0, 0, 0);
        for (let i = 0; i <= degree; i++) {
            const basisValue = this.basisFunctionDerivative(span - degree + i, degree, u, knots, order);
            vector = vector.add(ders[0][i].scale(basisValue * weights[span - degree + i]));
        }
        let weightSum = 0;
        for (let i = 0; i <= degree; i++) {
            weightSum += this.basisFunctionDerivative(span - degree + i, degree, u, knots, order) * weights[span - degree + i];
        }
        if (weightSum === 0) {
            return new Vector3D(0, 0, 0);
        }
        return vector.scale(1 / weightSum);
    }
    findSpan(n, degree, u, knots) {
        if (u >= knots[n + 1]) {
            return n;
        }
        if (u < knots[0]) {
            return 0;
        }
        let low = 0;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (u < knots[mid] || u >= knots[mid + 1]) {
            if (u < knots[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
}
export class Vertex {
    constructor(point) {
        this.point = point;
    }
    getPoint() {
        return this.point;
    }
    setPoint(point) {
        this.point = point;
    }
}
export class Edge {
    constructor(startVertex, endVertex, curve) {
        this.startVertex = startVertex;
        this.endVertex = endVertex;
        this.curve = curve;
    }
    getStartVertex() {
        return this.startVertex;
    }
    getEndVertex() {
        return this.endVertex;
    }
    getCurve() {
        return this.curve;
    }
    setCurve(curve) {
        this.curve = curve;
    }
    setStartVertex(startVertex) {
        this.startVertex = startVertex;
    }
    setEndVertex(endVertex) {
        this.endVertex = endVertex;
    }
    reverse() {
        const temp = this.startVertex;
        this.startVertex = this.endVertex;
        this.endVertex = temp;
        if (this.curve) {
            this.curve = this.curve.reverse();
        }
    }
    length() {
        if (!this.curve) {
            return this.startVertex.getPoint().distanceTo(this.endVertex.getPoint());
        }
        return this.curve.length();
    }
}
export class Face {
    constructor(surface, outerLoop, innerLoops) {
        this.surface = surface;
        this.outerLoop = outerLoop;
        this.innerLoops = innerLoops || [];
    }
    getSurface() {
        return this.surface;
    }
    getOuterLoop() {
        return this.outerLoop;
    }
    getInnerLoops() {
        return this.innerLoops;
    }
    setSurface(surface) {
        this.surface = surface;
    }
    addInnerLoop(loop) {
        this.innerLoops.push(loop);
    }
    removeInnerLoop(loop) {
        this.innerLoops = this.innerLoops.filter(innerLoop => innerLoop !== loop);
    }
}
export class Loop {
    constructor(edges) {
        this.edges = edges || [];
    }
    getEdges() {
        return this.edges;
    }
    addEdge(edge) {
        this.edges.push(edge);
    }
    removeEdge(edge) {
        this.edges = this.edges.filter(e => e !== edge);
    }
    isClosed() {
        if (this.edges.length < 2) {
            return false;
        }
        const firstVertex = this.edges[0].getStartVertex();
        const lastVertex = this.edges[this.edges.length - 1].getEndVertex();
        return firstVertex === lastVertex;
    }
}
export class Shell {
    constructor(faces) {
        this.faces = faces || [];
    }
    getFaces() {
        return this.faces;
    }
    addFace(face) {
        this.faces.push(face);
    }
    removeFace(face) {
        this.faces = this.faces.filter(f => f !== face);
    }
}
export class Solid {
    constructor(shells) {
        this.shells = shells || [];
    }
    getShells() {
        return this.shells;
    }
    addShell(shell) {
        this.shells.push(shell);
    }
    removeShell(shell) {
        this.shells = this.shells.filter(s => s !== shell);
    }
    volume() {
        let totalVolume = 0;
        for (const shell of this.shells) {
            for (const face of shell.getFaces()) {
                const surface = face.getSurface();
                if (!surface)
                    continue;
                const uMin = surface.getKnotsU()[0];
                const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
                const vMin = surface.getKnotsV()[0];
                const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
                const numSamplesU = 20;
                const numSamplesV = 20;
                const deltaU = (uMax - uMin) / numSamplesU;
                const deltaV = (vMax - vMin) / numSamplesV;
                for (let i = 0; i < numSamplesU; i++) {
                    for (let j = 0; j < numSamplesV; j++) {
                        const u = uMin + i * deltaU;
                        const v = vMin + j * deltaV;
                        const nextU = u + deltaU;
                        const nextV = v + deltaV;
                        const p0 = surface.evaluate(u, v);
                        const p1 = surface.evaluate(nextU, v);
                        const p2 = surface.evaluate(nextU, nextV);
                        const p3 = surface.evaluate(u, nextV);
                        if (p0 && p1 && p2 && p3) {
                            const v1 = new Vector3D(p1.getX() - p0.getX(), p1.getY() - p0.getY(), p1.getZ() - p0.getZ());
                            const v2 = new Vector3D(p2.getX() - p0.getX(), p2.getY() - p0.getY(), p2.getZ() - p0.getZ());
                            const v3 = new Vector3D(p3.getX() - p0.getX(), p3.getY() - p0.getY(), p3.getZ() - p0.getZ());
                            const tetVolume = 1 / 6 * Math.abs(v1.dot(v2.cross(v3)));
                            totalVolume += tetVolume;
                        }
                    }
                }
            }
        }
        return totalVolume;
    }
    surfaceArea() {
        let totalArea = 0;
        for (const shell of this.shells) {
            for (const face of shell.getFaces()) {
                const surface = face.getSurface();
                if (!surface)
                    continue;
                const uMin = surface.getKnotsU()[0];
                const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
                const vMin = surface.getKnotsV()[0];
                const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
                const numSamplesU = 20;
                const numSamplesV = 20;
                const deltaU = (uMax - uMin) / numSamplesU;
                const deltaV = (vMax - vMin) / numSamplesV;
                for (let i = 0; i < numSamplesU; i++) {
                    for (let j = 0; j < numSamplesV; j++) {
                        const u = uMin + i * deltaU;
                        const v = vMin + j * deltaV;
                        const nextU = u + deltaU;
                        const nextV = v + deltaV;
                        const p0 = surface.evaluate(u, v);
                        const p1 = surface.evaluate(nextU, v);
                        const p2 = surface.evaluate(u, nextV);
                        if (p0 && p1 && p2) {
                            const v1 = new Vector3D(p1.getX() - p0.getX(), p1.getY() - p0.getY(), p1.getZ() - p0.getZ());
                            const v2 = new Vector3D(p2.getX() - p0.getX(), p2.getY() - p0.getY(), p2.getZ() - p0.getZ());
                            const area = v1.cross(v2).magnitude() * 0.5;
                            totalArea += area;
                        }
                    }
                }
            }
        }
        return totalArea;
    }
    getFaces() {
        let faces = [];
        for (const shell of this.shells) {
            faces.push(...shell.getFaces());
        }
        return faces;
    }
    addFace(face) {
        for (const shell of this.shells) {
            shell.addFace(face);
        }
    }
    removeFace(face) {
        for (const shell of this.shells) {
            shell.removeFace(face);
        }
    }
    copy() {
        const newSolid = new Solid();
        for (const shell of this.shells) {
            const newShell = new Shell();
            for (const face of shell.getFaces()) {
                const newFace = new Face(face.surface, // Assuming surface is immutable or can be shared
                this._copyLoop(face.outerLoop), face.innerLoops.map(loop => this._copyLoop(loop)));
                newShell.addFace(newFace);
            }
            newSolid.addShell(newShell);
        }
        return newSolid;
    }
    _copyLoop(loop) {
        const newLoop = new Loop();
        for (const edge of loop.edges) {
            const newEdge = new Edge(new Vertex(edge.startVertex.point), new Vertex(edge.endVertex.point), edge.curve);
            newLoop.addEdge(newEdge);
        }
        return newLoop;
    }
}
export class Intersection {
    constructor(curve, surface1, surface2) {
        this.curve = curve;
        this.surface1 = surface1;
        this.surface2 = surface2;
    }
    getCurve() {
        return this.curve;
    }
    getSurface1() {
        return this.surface1;
    }
    getSurface2() {
        return this.surface2;
    }
}
export class SurfaceTrimmer {
    constructor() {
    }
    trimSurface(surface, intersection, keepSide) {
        if (!surface || !intersection || !keepSide) {
            return null;
        }
        const curve = intersection.getCurve();
        if (!curve) {
            return surface;
        }
        const points = [];
        const numSamples = 50;
        const start = curve.getKnots()[0];
        const end = curve.getKnots()[curve.getKnots().length - 1];
        const step = (end - start) / numSamples;
        for (let i = 0; i <= numSamples; i++) {
            const u = start + i * step;
            const point = curve.evaluate(u);
            if (point) {
                points.push(point);
            }
        }
        if (points.length < 2) {
            return surface;
        }
        const uvs = [];
        for (const point of points) {
            const uv = this._findUVOnSurface(surface, point);
            if (uv) {
                uvs.push(uv);
            }
        }
        if (uvs.length < 2)
            return surface;
        return this._trimSurfaceWithUVs(surface, uvs, keepSide);
    }
    splitSurface(surface, intersection) {
        if (!surface || !intersection) {
            return null;
        }
        const curve = intersection.getCurve();
        if (!curve) {
            return [surface];
        }
        const points = [];
        const numSamples = 50;
        const start = curve.getKnots()[0];
        const end = curve.getKnots()[curve.getKnots().length - 1];
        const step = (end - start) / numSamples;
        for (let i = 0; i <= numSamples; i++) {
            const u = start + i * step;
            const point = curve.evaluate(u);
            if (point) {
                points.push(point);
            }
        }
        if (points.length < 2) {
            return [surface];
        }
        const uvs = [];
        for (const point of points) {
            const uv = this._findUVOnSurface(surface, point);
            if (uv) {
                uvs.push(uv);
            }
        }
        if (uvs.length < 2) {
            return [surface];
        }
        return this._splitSurfaceWithUVs(surface, uvs);
    }
    isPointInside(surface, point) {
        if (!surface || !point) {
            return false;
        }
        const uv = this._findUVOnSurface(surface, point);
        if (!uv) {
            return false;
        }
        const uvs = [];
        const numSamples = 50;
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        const uStep = (uMax - uMin) / numSamples;
        const vStep = (vMax - vMin) / numSamples;
        for (let i = 0; i <= numSamples; i++) {
            for (let j = 0; j <= numSamples; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const p = surface.evaluate(u, v);
                if (p) {
                    const dist = p.distanceTo(point);
                    if (dist < 0.001) {
                        uvs.push({
                            u: u,
                            v: v
                        });
                    }
                }
            }
        }
        if (uvs.length === 0) {
            return false;
        }
        const testUv = uvs[0];
        if (!testUv)
            return false;
        return this._isPointInsideTrim(uvs, testUv, true);
    }
    _findUVOnSurface(surface, point) {
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        let u = (uMin + uMax) / 2;
        let v = (vMin + vMax) / 2;
        const tolerance = 0.001;
        const maxIterations = 20;
        for (let i = 0; i < maxIterations; i++) {
            const surfacePoint = surface.evaluate(u, v);
            if (!surfacePoint) {
                return null;
            }
            const diff = new Vector3D(surfacePoint.getX() - point.getX(), surfacePoint.getY() - point.getY(), surfacePoint.getZ() - point.getZ());
            if (diff.magnitude() < tolerance) {
                return {
                    u: u,
                    v: v
                };
            }
            const jacobian = this._calculateJacobian(surface, u, v);
            if (!jacobian) {
                return null;
            }
            const invJacobian = this._inverse2x2(jacobian);
            if (!invJacobian) {
                return null;
            }
            const deltaU = invJacobian[0][0] * diff.getX() + invJacobian[0][1] * diff.getY() + invJacobian[0][2] * diff.getZ();
            const deltaV = invJacobian[1][0] * diff.getX() + invJacobian[1][1] * diff.getY() + invJacobian[1][2] * diff.getZ();
            u -= deltaU;
            v -= deltaV;
            u = Math.max(uMin, Math.min(u, uMax));
            v = Math.max(vMin, Math.min(v, vMax));
        }
        return null;
    }
    _inverse2x2(matrix) {
        const det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        if (Math.abs(det) < 1e-8) {
            return null;
        }
        const invDet = 1 / det;
        return [
            [
                matrix[1][1] * invDet,
                -matrix[0][1] * invDet,
                0
            ],
            [
                -matrix[1][0] * invDet,
                matrix[0][0] * invDet,
                0
            ]
        ];
    }
    _trimSurfaceWithUVs(surface, uvs, keepSide) {
        if (!surface || !uvs || uvs.length < 2) {
            return surface;
        }
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        let newControlPoints = surface.getControlPoints().map(row => [...row]);
        let newWeights = surface.getWeights().map(row => [...row]);
        const numSamplesU = 20;
        const numSamplesV = 20;
        const uStep = (uMax - uMin) / numSamplesU;
        const vStep = (vMax - vMin) / numSamplesV;
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const point = surface.evaluate(u, v);
                if (!point)
                    continue;
                const isInside = this._isPointInsideTrim(uvs, {
                    u: u,
                    v: v
                }, keepSide);
                if (!isInside) {
                    for (let k = 0; k < newControlPoints.length; k++) {
                        for (let l = 0; l < newControlPoints[k].length; l++) {
                            const controlPoint = newControlPoints[k][l];
                            const weight = newWeights[k][l];
                            const cpUv = this._findUVOnSurface(surface, controlPoint);
                            if (cpUv) {
                                const dist = Math.sqrt((cpUv.u - u) * (cpUv.u - u) + (cpUv.v - v) * (cpUv.v - v));
                                if (dist < 0.1) {
                                    newControlPoints[k][l] = new Point3D(0, 0, 0);
                                    newWeights[k][l] = 0;
                                }
                            }
                        }
                    }
                }
            }
        }
        return new NurbsSurface(surface.getDegreeU(), surface.getDegreeV(), newControlPoints, surface.getKnotsU(), surface.getKnotsV(), newWeights);
    }
    _isPointInsideTrim(uvs, uv, keepSide) {
        if (!uvs || uvs.length < 2) {
            return false;
        }
        let windingNumber = 0;
        for (let i = 0; i < uvs.length - 1; i++) {
            const p1 = uvs[i];
            const p2 = uvs[i + 1];
            if (p1.v <= uv.v) {
                if (p2.v > uv.v) {
                    const vt = (uv.v - p1.v) / (p2.v - p1.v);
                    if (uv.u < p1.u + vt * (p2.u - p1.u)) {
                        windingNumber++;
                    }
                }
            } else {
                if (p2.v <= uv.v) {
                    const vt = (uv.v - p1.v) / (p2.v - p1.v);
                    if (uv.u < p1.u + vt * (p2.u - p1.u)) {
                        windingNumber--;
                    }
                }
            }
        }
        if (keepSide)
            return windingNumber !== 0;
        else
            return windingNumber === 0;
    }
    _splitSurfaceWithUVs(surface, uvs) {
        if (!surface || !uvs || uvs.length < 2) {
            return [surface];
        }
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        const numSamplesU = 20;
        const numSamplesV = 20;
        const uStep = (uMax - uMin) / numSamplesU;
        const vStep = (vMax - vMin) / numSamplesV;
        const surface1ControlPoints = surface.getControlPoints().map(row => [...row]);
        const surface1Weights = surface.getWeights().map(row => [...row]);
        const surface2ControlPoints = surface.getControlPoints().map(row => [...row]);
        const surface2Weights = surface.getWeights().map(row => [...row]);
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const point = surface.evaluate(u, v);
                if (!point)
                    continue;
                const isInside = this._isPointInsideTrim(uvs, {
                    u: u,
                    v: v
                }, true);
                if (isInside) {
                    for (let k = 0; k < surface1ControlPoints.length; k++) {
                        for (let l = 0; l < surface1ControlPoints[k].length; l++) {
                            const controlPoint = surface1ControlPoints[k][l];
                            const weight = surface1Weights[k][l];
                            const cpUv = this._findUVOnSurface(surface, controlPoint);
                            if (cpUv) {
                                const dist = Math.sqrt((cpUv.u - u) * (cpUv.u - u) + (cpUv.v - v) * (cpUv.v - v));
                                if (dist < 0.1) {
                                    surface2ControlPoints[k][l] = new Point3D(0, 0, 0);
                                    surface2Weights[k][l] = 0;
                                }
                            }
                        }
                    }
                } else {
                    for (let k = 0; k < surface2ControlPoints.length; k++) {
                        for (let l = 0; l < surface2ControlPoints[k].length; l++) {
                            const controlPoint = surface2ControlPoints[k][l];
                            const weight = surface2Weights[k][l];
                            const cpUv = this._findUVOnSurface(surface, controlPoint);
                            if (cpUv) {
                                const dist = Math.sqrt((cpUv.u - u) * (cpUv.u - u) + (cpUv.v - v) * (cpUv.v - v));
                                if (dist < 0.1) {
                                    surface1ControlPoints[k][l] = new Point3D(0, 0, 0);
                                    surface1Weights[k][l] = 0;
                                }
                            }
                        }
                    }
                }
            }
        }
        const trimmedSurface1 = new NurbsSurface(surface.getDegreeU(), surface.getDegreeV(), surface1ControlPoints, surface.getKnotsU(), surface.getKnotsV(), surface1Weights);
        const trimmedSurface2 = new NurbsSurface(surface.getDegreeU(), surface.getDegreeV(), surface2ControlPoints, surface.getKnotsU(), surface.getKnotsV(), surface2Weights);
        return [
            trimmedSurface1,
            trimmedSurface2
        ];
    }
}
export class BooleanOperation {
    constructor() {
        this.vertices = [];
    }
    union(solid1, solid2) {
        if (!solid1 || !solid2) {
            return null;
        }
        const solid1Copy = solid1.copy();
        const solid2Copy = solid2.copy();
        let allFaces1 = solid1Copy.getFaces();
        let allFaces2 = solid2Copy.getFaces();
        let intersections = [];
        //check for self intersections in each solid
        allFaces1 = allFaces1.reduce((acc, face) => {
            const surface = face.getSurface();
            if (this._hasSelfIntersections(surface)) {
                const splitSurfaces = this._splitSelfIntersectingSurface(surface);
                if (splitSurfaces && splitSurfaces.length > 0) {
                    splitSurfaces.forEach(splitSurface => {
                        acc.push(new Face(splitSurface, face.getOuterLoop(), face.getInnerLoops()));
                    });
                }
            } else {
                acc.push(face);
            }
            return acc;
        }, []);
        allFaces2 = allFaces2.reduce((acc, face) => {
            const surface = face.getSurface();
            if (this._hasSelfIntersections(surface)) {
                const splitSurfaces = this._splitSelfIntersectingSurface(surface);
                if (splitSurfaces && splitSurfaces.length > 0) {
                    splitSurfaces.forEach(splitSurface => {
                        acc.push(new Face(splitSurface, face.getOuterLoop(), face.getInnerLoops()));
                    });
                }
            } else {
                acc.push(face);
            }
            return acc;
        }, []);
        for (const face1 of allFaces1) {
            for (const face2 of allFaces2) {
                const surface1 = face1.getSurface();
                const surface2 = face2.getSurface();
                if (!surface1 || !surface2) {
                    continue;
                }
                if (this._isCoplanar(surface1, surface2)) {
                    if (this._isCoincident(surface1, surface2)) {
                        const trimmedSurface = this._handleCoplanarFaces(face1, face2, true);
                        if (trimmedSurface) {
                            const newFace = new Face(trimmedSurface, face1.getOuterLoop(), face1.getInnerLoops());
                            intersections.push(newFace);
                        }
                    } else {
                        const surfaceIntersections = this.findSurfaceIntersections(surface1, surface2);
                        if (surfaceIntersections && surfaceIntersections.length > 0) {
                            for (const intersection of surfaceIntersections) {
                                if (!this._isCurveOnBoundary(intersection.getCurve(), surface1) && !this._isCurveOnBoundary(intersection.getCurve(), surface2)) {
                                    intersections.push(intersection);
                                }
                            }
                        }
                    }
                    continue;
                }
                const surfaceIntersections = this.findSurfaceIntersections(surface1, surface2);
                if (surfaceIntersections && surfaceIntersections.length > 0) {
                    for (const intersection of surfaceIntersections) {
                        if (!this._isCurveOnBoundary(intersection.getCurve(), surface1) && !this._isCurveOnBoundary(intersection.getCurve(), surface2)) {
                            intersections.push(intersection);
                        }
                    }
                }
            }
        }
        // Handle cases where there are no intersections
        if (intersections.length === 0) {
            if (this._isSolidInside(solid1Copy, solid2Copy)) {
                return solid2Copy;
            } else if (this._isSolidInside(solid2Copy, solid1Copy)) {
                return solid1Copy;
            } else {
                const newSolid = new Solid();
                const newShell1 = new Shell(allFaces1);
                const newShell2 = new Shell(allFaces2);
                newSolid.addShell(newShell1);
                newSolid.addShell(newShell2);
                return newSolid;
            }
        }
        let allFaces = [
            ...allFaces1,
            ...allFaces2
        ];
        const trimmer = new SurfaceTrimmer();
        // Split surfaces based on intersections
        for (const intersection of intersections) {
            if (intersection instanceof Face)
                continue;
            const surface1 = intersection.getSurface1();
            const surface2 = intersection.getSurface2();
            allFaces = allFaces.reduce((acc, face) => {
                const surface = face.getSurface();
                if (surface === surface1 || surface === surface2) {
                    if (this._isCurveOnBoundary(intersection.getCurve(), surface)) {
                        acc.push(face);
                        return acc;
                    }
                    const splitSurfaces = trimmer.splitSurface(surface, intersection);
                    if (splitSurfaces && splitSurfaces.length > 0) {
                        splitSurfaces.forEach(splitSurface => {
                            if (!this._hasSelfIntersections(splitSurface)) {
                                acc.push(new Face(splitSurface, face.getOuterLoop(), face.getInnerLoops()));
                            }
                        });
                    } else {
                        acc.push(face);
                    }
                } else {
                    acc.push(face);
                }
                return acc;
            }, []);
        }
        //remove any faces that have a surface with zero area
        allFaces = allFaces.filter(face => {
            const surface = face.getSurface();
            if (!surface)
                return false;
            const uMin = surface.getKnotsU()[0];
            const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
            const vMin = surface.getKnotsV()[0];
            const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
            const area = (uMax - uMin) * (vMax - vMin);
            if (area <= 0) {
                return false;
            }
            return true;
        });
        // Filter out tiny faces
        allFaces = allFaces.filter(face => this._isFaceLargeEnough(face));
        const newSolid = new Solid();
        const newShell = new Shell(allFaces);
        newSolid.addShell(newShell);
        return newSolid;
    }
    intersect(solid1, solid2) {
        if (!solid1 || !solid2) {
            return null;
        }
        const solid1Copy = solid1.copy();
        const solid2Copy = solid2.copy();
        const allFaces1 = solid1Copy.getFaces();
        const allFaces2 = solid2Copy.getFaces();
        const intersections = [];
        // Check if one solid is fully inside the other
        if (this._isSolidInside(solid1Copy, solid2Copy)) {
            return solid1Copy;
        }
        if (this._isSolidInside(solid2Copy, solid1Copy)) {
            return solid2Copy;
        }
        for (const face1 of allFaces1) {
            for (const face2 of allFaces2) {
                const surface1 = face1.getSurface();
                const surface2 = face2.getSurface();
                if (!surface1 || !surface2) {
                    continue;
                }
                const surfaceIntersections = this.findSurfaceIntersections(surface1, surface2);
                if (this._isCoplanar(surface1, surface2)) {
                    if (surfaceIntersections && surfaceIntersections.length > 0) {
                        for (const intersection of surfaceIntersections) {
                            if (!this._isCurveOnBoundary(intersection.getCurve(), surface1) && !this._isCurveOnBoundary(intersection.getCurve(), surface2)) {
                                intersections.push(intersection);
                            }
                        }
                    } else {
                        const trimmedSurface = this._handleCoplanarFaces(face1, face2, true);
                        if (trimmedSurface) {
                            const newFace = new Face(trimmedSurface, face1.getOuterLoop(), face1.getInnerLoops());
                            intersections.push(newFace);
                        }
                    }
                    continue;
                }
                if (surfaceIntersections && surfaceIntersections.length > 0) {
                    for (const intersection of surfaceIntersections) {
                        if (!this._isCurveOnBoundary(intersection.getCurve(), surface1) && !this._isCurveOnBoundary(intersection.getCurve(), surface2)) {
                            intersections.push(intersection);
                        }
                    }
                }
            }
        }
        const allFaces = [
            ...allFaces1,
            ...allFaces2
        ];
        let trimmedFaces = [];
        const trimmer = new SurfaceTrimmer();
        for (const face of allFaces) {
            let currentSurface = face.getSurface();
            if (!currentSurface)
                continue;
            let splitFaces = [face];
            for (const intersection of intersections) {
                if (intersection instanceof Face) {
                    if (intersection.getSurface() === currentSurface) {
                        splitFaces = [intersection];
                        currentSurface = null;
                        break;
                    }
                    continue;
                }
                const surface1 = intersection.getSurface1();
                const surface2 = intersection.getSurface2();
                if (surface1 === currentSurface) {
                    let newSplitFaces = [];
                    for (const splitFace of splitFaces) {
                        const splitSurfaces = trimmer.splitSurface(splitFace.getSurface(), intersection);
                        if (splitSurfaces && splitSurfaces.length > 0) {
                            splitSurfaces.forEach(splitSurface => {
                                newSplitFaces.push(new Face(splitSurface, splitFace.getOuterLoop(), splitFace.getInnerLoops()));
                            });
                        } else {
                            newSplitFaces.push(splitFace);
                        }
                    }
                    splitFaces = newSplitFaces;
                    if (splitFaces.length > 0)
                        currentSurface = splitFaces[0].getSurface();
                    else
                        currentSurface = null;
                } else if (surface2 === currentSurface) {
                    let newSplitFaces = [];
                    for (const splitFace of splitFaces) {
                        const splitSurfaces = trimmer.splitSurface(splitFace.getSurface(), intersection);
                        if (splitSurfaces && splitSurfaces.length > 0) {
                            splitSurfaces.forEach(splitSurface => {
                                newSplitFaces.push(new Face(splitSurface, splitFace.getOuterLoop(), splitFace.getInnerLoops()));
                            });
                        } else {
                            newSplitFaces.push(splitFace);
                        }
                    }
                    splitFaces = newSplitFaces;
                    if (splitFaces.length > 0)
                        currentSurface = splitFaces[0].getSurface();
                    else
                        currentSurface = null;
                }
            }
            if (splitFaces && splitFaces.length > 0) {
                trimmedFaces.push(...splitFaces);
            }
        }
        // Filter out tiny faces
        trimmedFaces = trimmedFaces.filter(face => this._isFaceLargeEnough(face));
        const newSolid = new Solid();
        const newShell = new Shell(trimmedFaces);
        newSolid.addShell(newShell);
        return newSolid;
    }
    subtract(solid1, solid2) {
        if (!solid1 || !solid2) {
            return null;
        }
        const solid1Copy = solid1.copy();
        const solid2Copy = solid2.copy();
        const allFaces1 = solid1Copy.getFaces();
        const allFaces2 = solid2Copy.getFaces();
        const intersections = [];
        // Check if solid2 is fully inside solid1
        if (this._isSolidInside(solid2Copy, solid1Copy)) {
            const newSolid = new Solid();
            const newShell = new Shell(allFaces1);
            newSolid.addShell(newShell);
            return newSolid;
        }
        if (this._isSolidInside(solid1Copy, solid2Copy)) {
            return new Solid();
        }
        for (const face1 of allFaces1) {
            for (const face2 of allFaces2) {
                const surface1 = face1.getSurface();
                const surface2 = face2.getSurface();
                if (!surface1 || !surface2) {
                    continue;
                }
                const surfaceIntersections = this.findSurfaceIntersections(surface1, surface2);
                if (this._isCoplanar(surface1, surface2)) {
                    if (surfaceIntersections && surfaceIntersections.length > 0) {
                        for (const intersection of surfaceIntersections) {
                            if (!this._isCurveOnBoundary(intersection.getCurve(), surface1) && !this._isCurveOnBoundary(intersection.getCurve(), surface2)) {
                                intersections.push(intersection);
                            }
                        }
                    } else {
                        const trimmedSurface = this._handleCoplanarFaces(face1, face2, false);
                        if (trimmedSurface) {
                            const newFace = new Face(trimmedSurface, face1.getOuterLoop(), face1.getInnerLoops());
                            intersections.push(newFace);
                        }
                    }
                    continue;
                }
                if (surfaceIntersections && surfaceIntersections.length > 0) {
                    for (const intersection of surfaceIntersections) {
                        if (!this._isCurveOnBoundary(intersection.getCurve(), surface1) && !this._isCurveOnBoundary(intersection.getCurve(), surface2)) {
                            intersections.push(intersection);
                        }
                    }
                }
            }
        }
        const allFaces = [
            ...allFaces1,
            ...allFaces2
        ];
        let trimmedFaces = [];
        const trimmer = new SurfaceTrimmer();
        for (const face of allFaces) {
            let currentSurface = face.getSurface();
            if (!currentSurface)
                continue;
            let splitFaces = [face];
            for (const intersection of intersections) {
                if (intersection instanceof Face) {
                    if (intersection.getSurface() === currentSurface) {
                        splitFaces = [intersection];
                        currentSurface = null;
                        break;
                    }
                    continue;
                }
                const surface1 = intersection.getSurface1();
                const surface2 = intersection.getSurface2();
                if (allFaces1.includes(face) && surface1 === currentSurface) {
                    let newSplitFaces = [];
                    for (const splitFace of splitFaces) {
                        const splitSurfaces = trimmer.splitSurface(splitFace.getSurface(), intersection);
                        if (splitSurfaces && splitSurfaces.length > 0) {
                            splitSurfaces.forEach(splitSurface => {
                                newSplitFaces.push(new Face(splitSurface, splitFace.getOuterLoop(), splitFace.getInnerLoops()));
                            });
                        } else {
                            newSplitFaces.push(splitFace);
                        }
                    }
                    splitFaces = newSplitFaces;
                    if (splitFaces.length > 0)
                        currentSurface = splitFaces[0].getSurface();
                    else
                        currentSurface = null;
                } else if (allFaces2.includes(face) && surface2 === currentSurface) {
                    let newSplitFaces = [];
                    for (const splitFace of splitFaces) {
                        const splitSurfaces = trimmer.splitSurface(splitFace.getSurface(), intersection);
                        if (splitSurfaces && splitSurfaces.length > 0) {
                            splitSurfaces.forEach(splitSurface => {
                                newSplitFaces.push(new Face(splitSurface, splitFace.getOuterLoop(), splitFace.getInnerLoops()));
                            });
                        } else {
                            newSplitFaces.push(splitFace);
                        }
                    }
                    splitFaces = newSplitFaces;
                    if (splitFaces.length > 0)
                        currentSurface = splitFaces[0].getSurface();
                    else
                        currentSurface = null;
                }
            }
            if (splitFaces && splitFaces.length > 0) {
                trimmedFaces.push(...splitFaces);
            }
        }
        // Filter out tiny faces
        trimmedFaces = trimmedFaces.filter(face => this._isFaceLargeEnough(face));
        const newSolid = new Solid();
        const newShell = new Shell(trimmedFaces);
        newSolid.addShell(newShell);
        return newSolid;
    }
    findSurfaceIntersections(surface1, surface2) {
        if (!surface1 || !surface2) {
            return null;
        }
        if (this._isCoplanar(surface1, surface2)) {
            if (this._isCoincident(surface1, surface2)) {
                return [];
            }
        }
        const initialCurves = this._approximateIntersections(surface1, surface2);
        if (!initialCurves || initialCurves.length === 0) {
            return [];
        }
        const refinedCurves = [];
        const tolerance = 0.001;
        const maxIterations = 20;
        for (const initialCurve of initialCurves) {
            const refinedPoints = [];
            for (const initialPoint of initialCurve) {
                const refinedPoint = this._refineIntersectionPoint(surface1, surface2, initialPoint, tolerance, maxIterations);
                if (refinedPoint) {
                    refinedPoints.push(refinedPoint);
                }
            }
            if (refinedPoints.length > 1) {
                const controlPoints = refinedPoints.map(point => surface1.evaluate(point.u1, point.v1));
                if (controlPoints && controlPoints.length > 1 && !controlPoints.includes(null)) {
                    const degree = 3;
                    const knots = this._createUniformKnotVector(controlPoints.length, degree);
                    const weights = Array(controlPoints.length).fill(1);
                    const nurbsCurve = new NurbsCurve(degree, controlPoints, knots, weights);
                    // Check if the refined curve is within the bounds of both surfaces.
                    if (this._isCurveWithinBounds(nurbsCurve, surface1) && this._isCurveWithinBounds(nurbsCurve, surface2)) {
                        refinedCurves.push(nurbsCurve);
                    }
                }
            }
        }
        const filteredCurves = refinedCurves.filter(curve => curve.length() > 0.01);
        const intersections = filteredCurves.map(curve => new Intersection(curve, surface1, surface2));
        const validIntersections = [];
        if (this._isCoplanar(surface1, surface2)) {
            return intersections;
        }
        for (const intersection of intersections) {
            if (!this._isCurveOnBoundary(intersection.getCurve(), surface1) && !this._isCurveOnBoundary(intersection.getCurve(), surface2)) {
                validIntersections.push(intersection);
            }
        }
        return validIntersections;
    }
    _findSurfaceIntersections(surface1, surface2) {
        if (!surface1 || !surface2) {
            return null;
        }
        const initialCurves = this._approximateIntersections(surface1, surface2);
        if (!initialCurves || initialCurves.length === 0) {
            return [];
        }
        const refinedCurves = this._refineIntersections(surface1, surface2, initialCurves);
        if (!refinedCurves || refinedCurves.length === 0) {
            return [];
        }
        const intersections = refinedCurves.map(curve => new Intersection(curve, surface1, surface2));
        return intersections;
    }
    _approximateIntersections(surface1, surface2) {
        if (!surface1 || !surface2) {
            return null;
        }
        const uMin1 = surface1.getKnotsU()[0];
        const uMax1 = surface1.getKnotsU()[surface1.getKnotsU().length - 1];
        const vMin1 = surface1.getKnotsV()[0];
        const vMax1 = surface1.getKnotsV()[surface1.getKnotsV().length - 1];
        const uMin2 = surface2.getKnotsU()[0];
        const uMax2 = surface2.getKnotsU()[surface2.getKnotsU().length - 1];
        const vMin2 = surface2.getKnotsV()[0];
        const vMax2 = surface2.getKnotsV()[surface2.getKnotsV().length - 1];
        const initialCurves = [];
        this._recursiveSubdivision(surface1, surface2, uMin1, uMax1, vMin1, vMax1, uMin2, uMax2, vMin2, vMax2, initialCurves, 0);
        return initialCurves;
    }
    _refineIntersections(surface1, surface2, initialCurves) {
        if (!surface1 || !surface2 || !initialCurves || initialCurves.length === 0) {
            return null;
        }
        const refinedCurves = [];
        const tolerance = 0.001;
        const maxIterations = 20;
        for (const initialCurve of initialCurves) {
            const refinedPoints = [];
            for (const initialPoint of initialCurve) {
                const refinedPoint = this._refineIntersectionPoint(surface1, surface2, initialPoint, tolerance, maxIterations);
                if (refinedPoint) {
                    refinedPoints.push(refinedPoint);
                }
            }
            if (refinedPoints.length > 1) {
                const controlPoints = refinedPoints.map(point => surface1.evaluate(point.u1, point.v1));
                if (controlPoints && controlPoints.length > 1) {
                    const degree = 3;
                    const knots = this._createUniformKnotVector(controlPoints.length, degree);
                    const weights = Array(controlPoints.length).fill(1);
                    const nurbsCurve = new NurbsCurve(degree, controlPoints, knots, weights);
                    refinedCurves.push(nurbsCurve);
                }
            }
        }
        return refinedCurves;
    }
    _refineIntersectionPoint(surface1, surface2, initialPoint, tolerance, maxIterations) {
        let u1 = initialPoint.u1;
        let v1 = initialPoint.v1;
        let u2 = initialPoint.u2;
        let v2 = initialPoint.v2;
        const refinedPoint = this._levenbergMarquardtRefinement(surface1, surface2, {
            u1,
            v1,
            u2,
            v2
        }, tolerance, maxIterations);
        if (refinedPoint)
            return refinedPoint;
        return null;
    }
    _calculateJacobian(surface, u, v) {
        const du = surface.evaluatePartialDerivative(u, v, 1, 0);
        const dv = surface.evaluatePartialDerivative(u, v, 0, 1);
        if (!du || !dv) {
            return null;
        }
        return [
            [
                du.getX(),
                du.getY(),
                du.getZ()
            ],
            [
                dv.getX(),
                dv.getY(),
                dv.getZ()
            ]
        ];
    }
    _createUniformKnotVector(numPoints, degree) {
        if (numPoints <= degree + 1) {
            return null;
        }
        const knots = [];
        for (let i = 0; i <= degree; i++) {
            knots.push(0);
        }
        for (let i = 1; i <= numPoints - degree - 1; i++) {
            knots.push(i);
        }
        for (let i = 0; i <= degree; i++) {
            knots.push(numPoints - degree - 1);
        }
        return knots;
    }
    _createEdgesFromIntersection(intersection, face1, face2) {
        if (!intersection || !face1 || !face2) {
            return null;
        }
        let curve = intersection.getCurve();
        if (!curve || curve.getControlPoints().length < 2) {
            return [];
        }
        const edges = [];
        // Split the curve at self-intersections
        const numSamples = 50;
        const start = curve.getKnots()[0];
        const end = curve.getKnots()[curve.getKnots().length - 1];
        const step = (end - start) / numSamples;
        const points = [];
        for (let i = 0; i <= numSamples; i++) {
            const u = start + i * step;
            const point = curve.evaluate(u);
            if (point) {
                points.push({
                    point: point,
                    u: u
                });
            }
        }
        if (points.length < 2) {
            return [];
        }
        let splitCurves = [curve];
        let newSplitCurves = [];
        for (let i = 0; i < points.length; i++) {
            for (let j = i + 1; j < points.length; j++) {
                const p1 = points[i].point;
                const p2 = points[j].point;
                if (p1.distanceTo(p2) < 0.001) {
                    if (Math.abs(points[i].u - points[j].u) > 0.001) {
                        newSplitCurves = [];
                        for (const splitCurve of splitCurves) {
                            const splitResult = splitCurve.splitCurve(points[j].u);
                            if (splitResult && splitResult.length === 2) {
                                newSplitCurves.push(splitResult[0]);
                                newSplitCurves.push(splitResult[1]);
                            } else {
                                newSplitCurves.push(splitCurve);
                            }
                        }
                        splitCurves = newSplitCurves;
                    }
                }
            }
        }
        for (const splitCurve of splitCurves) {
            const numSamples = 50;
            const start = splitCurve.getKnots()[0];
            const end = splitCurve.getKnots()[splitCurve.getKnots().length - 1];
            const step = (end - start) / numSamples;
            const points = [];
            for (let i = 0; i <= numSamples; i++) {
                const u = start + i * step;
                const point = splitCurve.evaluate(u);
                if (point) {
                    points.push(point);
                }
            }
            if (points.length < 2) {
                continue;
            }
            for (let i = 0; i < points.length - 1; i++) {
                const startPoint = points[i];
                const endPoint = points[i + 1];
                const startVertex = this._findOrCreateVertex(startPoint);
                const endVertex = this._findOrCreateVertex(endPoint);
                if (startVertex === endVertex) {
                    continue;
                }
                const edgeCurve = new NurbsCurve(splitCurve.getDegree(), [
                    startPoint,
                    endPoint
                ], [
                    0,
                    1
                ], [
                    1,
                    1
                ]);
                const edge = new Edge(startVertex, endVertex, edgeCurve);
                if (edge.length() > 0.001) {
                    edges.push(edge);
                }
            }
        }
        return edges;
    }
    _findOrCreateVertex(point) {
        if (!this.vertices) {
            this.vertices = [];
        }
        for (const vertex of this.vertices) {
            if (vertex.getPoint().distanceTo(point) < 0.001) {
                return vertex;
            }
        }
        const newVertex = new Vertex(point);
        this.vertices.push(newVertex);
        return newVertex;
    }
    _createLoopsFromEdges(edges) {
        if (!edges || edges.length === 0) {
            return [];
        }
        const loops = [];
        let remainingEdges = [...edges];
        while (remainingEdges.length > 0) {
            const currentEdge = remainingEdges.shift();
            const currentLoop = new Loop([currentEdge]);
            let currentVertex = currentEdge.getEndVertex();
            let loopComplete = false;
            while (!loopComplete) {
                let nextEdgeIndex = -1;
                let minDistance = Infinity;
                for (let i = 0; i < remainingEdges.length; i++) {
                    const nextEdge = remainingEdges[i];
                    const startVertex = nextEdge.getStartVertex();
                    const distance = currentVertex.getPoint().distanceTo(startVertex.getPoint());
                    if (distance < minDistance) {
                        minDistance = distance;
                        nextEdgeIndex = i;
                    }
                }
                if (nextEdgeIndex !== -1 && minDistance < 0.001) {
                    const nextEdge = remainingEdges.splice(nextEdgeIndex, 1)[0];
                    let isDuplicate = false;
                    for (const edge of currentLoop.getEdges()) {
                        if (edge === nextEdge) {
                            isDuplicate = true;
                            break;
                        }
                    }
                    if (!isDuplicate) {
                        let edgeIsInLoop = false;
                        for (const edge of currentLoop.getEdges()) {
                            if (edge.getStartVertex().getPoint().distanceTo(nextEdge.getStartVertex().getPoint()) < 0.001 && edge.getEndVertex().getPoint().distanceTo(nextEdge.getEndVertex().getPoint()) < 0.001) {
                                edgeIsInLoop = true;
                                break;
                            }
                        }
                        if (!edgeIsInLoop) {
                            currentLoop.addEdge(nextEdge);
                            currentVertex = nextEdge.getEndVertex();
                        } else {
                            loopComplete = true;
                        }
                    } else {
                        loopComplete = true;
                    }
                } else {
                    let closingEdgeIndex = -1;
                    for (let i = 0; i < currentLoop.getEdges().length; i++) {
                        if (currentLoop.getEdges()[i].getStartVertex().getPoint().distanceTo(currentVertex.getPoint()) < 0.001) {
                            closingEdgeIndex = i;
                            break;
                        }
                    }
                    if (closingEdgeIndex !== -1) {
                        currentLoop.addEdge(currentLoop.getEdges()[closingEdgeIndex]);
                    }
                    loopComplete = true;
                }
                if (currentLoop.isClosed()) {
                    loopComplete = true;
                }
            }
            if (currentLoop.isClosed() && this._isLoopValid(currentLoop)) {
                loops.push(currentLoop);
            }
        }
        return loops;
    }
    _isCoincident(surface1, surface2) {
        if (!surface1 || !surface2) {
            return false;
        }
        const uMin1 = surface1.getKnotsU()[0];
        const uMax1 = surface1.getKnotsU()[surface1.getKnotsU().length - 1];
        const vMin1 = surface1.getKnotsV()[0];
        const vMax1 = surface1.getKnotsV()[surface1.getKnotsV().length - 1];
        const uMin2 = surface2.getKnotsU()[0];
        const uMax2 = surface2.getKnotsU()[surface2.getKnotsU().length - 1];
        const vMin2 = surface2.getKnotsV()[0];
        const vMax2 = surface2.getKnotsV()[surface2.getKnotsV().length - 1];
        //check for overlapping parameter domains
        if (uMax1 < uMin2 || uMin1 > uMax2 || vMax1 < vMin2 || vMin1 > vMax2) {
            return false;
        }
        const numSamplesU = 10;
        const numSamplesV = 10;
        const deltaU1 = (uMax1 - uMin1) / numSamplesU;
        const deltaV1 = (vMax1 - vMin1) / numSamplesV;
        const tolerance = 0.001;
        let hasOverlap = false;
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u1 = uMin1 + i * deltaU1;
                const v1 = vMin1 + j * deltaV1;
                const point1 = surface1.evaluate(u1, v1);
                if (!point1)
                    continue;
                const uv2 = this._findUVOnSurface(surface2, point1);
                if (uv2) {
                    const point2 = surface2.evaluate(uv2.u, uv2.v);
                    if (!point2)
                        continue;
                    const dist = point1.distanceTo(point2);
                    if (dist <= tolerance) {
                        hasOverlap = true;
                    }
                }
            }
        }
        return hasOverlap;
    }
    _isTangentIntersection(intersection, surface1, surface2) {
        if (!intersection || !surface1 || !surface2) {
            return false;
        }
        const curve = intersection.getCurve();
        if (!curve) {
            return false;
        }
        const numSamples = 10;
        const start = curve.getKnots()[0];
        const end = curve.getKnots()[curve.getKnots().length - 1];
        const step = (end - start) / numSamples;
        for (let i = 0; i <= numSamples; i++) {
            const u = start + i * step;
            const point = curve.evaluate(u);
            if (!point)
                continue;
            const uv1 = this._findUVOnSurface(surface1, point);
            const uv2 = this._findUVOnSurface(surface2, point);
            if (!uv1 || !uv2)
                continue;
            const normal1 = surface1.normal(uv1.u, uv1.v);
            const normal2 = surface2.normal(uv2.u, uv2.v);
            if (!normal1 || !normal2)
                continue;
            const dotProduct = normal1.dot(normal2);
            if (Math.abs(dotProduct) > 0.99) {
                return true;
            }
        }
        return false;
    }
    _recursiveSubdivision(surface1, surface2, uMin1, uMax1, vMin1, vMax1, uMin2, uMax2, vMin2, vMax2, initialCurves, depth) {
        const maxDepth = 10;
        // Increased max depth
        if (depth > maxDepth) {
            return;
        }
        const uMid1 = (uMin1 + uMax1) / 2;
        const vMid1 = (vMin1 + vMax1) / 2;
        const uMid2 = (uMin2 + uMax2) / 2;
        const vMid2 = (vMin2 + vMax2) / 2;
        const p1Min = surface1.evaluate(uMin1, vMin1);
        const p1Max = surface1.evaluate(uMax1, vMax1);
        const p1Mid = surface1.evaluate(uMid1, vMid1);
        const p2Min = surface2.evaluate(uMin2, vMin2);
        const p2Max = surface2.evaluate(uMax2, vMax2);
        const p2Mid = surface2.evaluate(uMid2, vMid2);
        if (!p1Min || !p1Max || !p1Mid || !p2Min || !p2Max || !p2Mid)
            return;
        const bbox1 = this._getBoundingBox(surface1, uMin1, uMax1, vMin1, vMax1);
        const bbox2 = this._getBoundingBox(surface2, uMin2, uMax2, vMin2, vMax2);
        if (!this._doBoundingBoxesOverlap(bbox1, bbox2)) {
            return;
        }
        const distMin = p1Min.distanceTo(p2Min);
        const distMax = p1Max.distanceTo(p2Max);
        const distMid = p1Mid.distanceTo(p2Mid);
        const threshold = 0.1;
        if (distMin < threshold && distMax < threshold && distMid < threshold) {
            initialCurves.push([{
                    u1: uMid1,
                    v1: vMid1,
                    u2: uMid2,
                    v2: vMid2
                }]);
            return;
        }
        // More aggressive subdivision
        const distThreshold = 0.5;
        if (distMin < distThreshold || distMax < distThreshold || distMid < distThreshold) {
            this._recursiveSubdivision(surface1, surface2, uMin1, uMid1, vMin1, vMid1, uMin2, uMid2, vMin2, vMid2, initialCurves, depth + 1);
            this._recursiveSubdivision(surface1, surface2, uMid1, uMax1, vMin1, vMid1, uMid2, uMax2, vMin2, vMid2, initialCurves, depth + 1);
            this._recursiveSubdivision(surface1, surface2, uMin1, uMid1, vMid1, vMax1, uMin2, uMid2, vMid2, vMax2, initialCurves, depth + 1);
            this._recursiveSubdivision(surface1, surface2, uMid1, uMax1, vMid1, vMax1, uMid2, uMax2, vMid2, vMax2, initialCurves, depth + 1);
        } else {
            this._recursiveSubdivision(surface1, surface2, uMin1, uMid1, vMin1, vMid1, uMin2, uMid2, vMin2, vMid2, initialCurves, depth + 1);
            this._recursiveSubdivision(surface1, surface2, uMid1, uMax1, vMin1, vMid1, uMid2, uMax2, vMin2, vMid2, initialCurves, depth + 1);
            this._recursiveSubdivision(surface1, surface2, uMin1, uMid1, vMid1, vMax1, uMin2, uMid2, vMid2, vMax2, initialCurves, depth + 1);
            this._recursiveSubdivision(surface1, surface2, uMid1, uMax1, vMid1, vMax1, uMid2, uMax2, vMid2, vMax2, initialCurves, depth + 1);
        }
    }
    _isFaceInside(surface1, surface2) {
        if (!surface1 || !surface2) {
            return false;
        }
        const numSamplesU = 5;
        const numSamplesV = 5;
        const uMin1 = surface1.getKnotsU()[0];
        const uMax1 = surface1.getKnotsU()[surface1.getKnotsU().length - 1];
        const vMin1 = surface1.getKnotsV()[0];
        const vMax1 = surface1.getKnotsV()[surface1.getKnotsV().length - 1];
        const uStep1 = (uMax1 - uMin1) / numSamplesU;
        const vStep1 = (vMax1 - vMin1) / numSamplesV;
        let allInside = true;
        const trimmer = new SurfaceTrimmer();
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u1 = uMin1 + i * uStep1;
                const v1 = vMin1 + j * vStep1;
                const point1 = surface1.evaluate(u1, v1);
                if (point1) {
                    if (!trimmer.isPointInside(surface2, point1)) {
                        allInside = false;
                    }
                }
            }
        }
        return allInside;
    }
    _hasSelfIntersections(surface) {
        if (!surface) {
            return false;
        }
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        const numSamplesU = 20;
        const numSamplesV = 20;
        const uStep = (uMax - uMin) / numSamplesU;
        const vStep = (vMax - vMin) / numSamplesV;
        const points = [];
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const point = surface.evaluate(u, v);
                if (point)
                    points.push({
                        u,
                        v,
                        point
                    });
            }
        }
        for (let i = 0; i < points.length; i++) {
            for (let j = i + 1; j < points.length; j++) {
                const point1 = points[i].point;
                const point2 = points[j].point;
                const dist = point1.distanceTo(point2);
                if (dist < 0.001) {
                    const u1 = points[i].u;
                    const v1 = points[i].v;
                    const u2 = points[j].u;
                    const v2 = points[j].v;
                    if (Math.abs(u1 - u2) > 0.001 || Math.abs(v1 - v2) > 0.001) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
    _isCurveOnBoundary(curve, surface) {
        if (!curve || !surface) {
            return false;
        }
        const numSamples = 50;
        const start = curve.getKnots()[0];
        const end = curve.getKnots()[curve.getKnots().length - 1];
        const step = (end - start) / numSamples;
        for (let i = 0; i <= numSamples; i++) {
            const u = start + i * step;
            const point = curve.evaluate(u);
            if (!point)
                continue;
            const uv = this._findUVOnSurface(surface, point);
            if (!uv)
                return false;
            const uMin = surface.getKnotsU()[0];
            const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
            const vMin = surface.getKnotsV()[0];
            const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
            const tolerance = 0.001;
            if (Math.abs(uv.u - uMin) < tolerance || Math.abs(uv.u - uMax) < tolerance || Math.abs(uv.v - vMin) < tolerance || Math.abs(uv.v - vMax) < tolerance) {
                continue;
            } else {
                return false;
            }
        }
        return true;
    }
    _levenbergMarquardtRefinement(surface1, surface2, initialPoint, tolerance, maxIterations) {
        let u1 = initialPoint.u1;
        let v1 = initialPoint.v1;
        let u2 = initialPoint.u2;
        let v2 = initialPoint.v2;
        let lambda = 0.01;
        const dampingFactor = 10;
        for (let i = 0; i < maxIterations; i++) {
            const p1 = surface1.evaluate(u1, v1);
            const p2 = surface2.evaluate(u2, v2);
            if (!p1 || !p2) {
                return null;
            }
            const diff = new Vector3D(p2.getX() - p1.getX(), p2.getY() - p1.getY(), p2.getZ() - p1.getZ());
            if (diff.magnitude() < tolerance) {
                return {
                    u1: u1,
                    v1: v1,
                    u2: u2,
                    v2: v2
                };
            }
            const jacobian = this._calculateJacobianForIntersection(surface1, surface2, u1, v1, u2, v2);
            if (!jacobian) {
                return null;
            }
            const JtJ = this._transposeMultiply(jacobian, jacobian);
            const JtDiff = this._transposeMultiply(jacobian, [
                [diff.getX()],
                [diff.getY()],
                [diff.getZ()]
            ]);
            if (!JtJ || !JtDiff) {
                return null;
            }
            const augmentedJtJ = this._augmentMatrix(JtJ, lambda);
            const delta = this._solveLinearSystem(augmentedJtJ, JtDiff);
            if (!delta) {
                lambda *= dampingFactor;
                continue;
            }
            if (isNaN(delta[0][0]) || isNaN(delta[1][0]) || isNaN(delta[2][0]) || isNaN(delta[3][0]))
                return null;
            const prevU1 = u1;
            const prevV1 = v1;
            const prevU2 = u2;
            const prevV2 = v2;
            u1 -= delta[0][0];
            v1 -= delta[1][0];
            u2 -= delta[2][0];
            v2 -= delta[3][0];
            const uMin1 = surface1.getKnotsU()[0];
            const uMax1 = surface1.getKnotsU()[surface1.getKnotsU().length - 1];
            const vMin1 = surface1.getKnotsV()[0];
            const vMax1 = surface1.getKnotsV()[surface1.getKnotsV().length - 1];
            const uMin2 = surface2.getKnotsU()[0];
            const uMax2 = surface2.getKnotsU()[surface2.getKnotsU().length - 1];
            const vMin2 = surface2.getKnotsV()[0];
            const vMax2 = surface2.getKnotsV()[surface2.getKnotsV().length - 1];
            u1 = Math.max(uMin1, Math.min(u1, uMax1));
            v1 = Math.max(vMin1, Math.min(v1, vMax1));
            u2 = Math.max(uMin2, Math.min(u2, uMax2));
            v2 = Math.max(vMin2, Math.min(v2, vMax2));
            if (isNaN(u1) || isNaN(v1) || isNaN(u2) || isNaN(v2)) {
                return null;
            }
            if (Math.abs(delta[0][0]) < 1e-10 && Math.abs(delta[1][0]) < 1e-10 && Math.abs(delta[2][0]) < 1e-10 && Math.abs(delta[3][0]) < 1e-10) {
                return null;
            }
            if (delta[0][0] === 0 && delta[1][0] === 0 && delta[2][0] === 0 && delta[3][0] === 0) {
                lambda *= 0.1;
                continue;
            }
            const stepSize = Math.sqrt((u1 - prevU1) * (u1 - prevU1) + (v1 - prevV1) * (v1 - prevV1) + (u2 - prevU2) * (u2 - prevU2) + (v2 - prevV2) * (v2 - prevV2));
            if (stepSize > 0.1) {
                u1 = prevU1 - delta[0][0] * 0.1 / stepSize;
                v1 = prevV1 - delta[1][0] * 0.1 / stepSize;
                u2 = prevU2 - delta[2][0] * 0.1 / stepSize;
                v2 = prevV2 - delta[3][0] * 0.1 / stepSize;
            }
            lambda *= 0.1;
            if (lambda < 1e-10) {
                lambda = 1e-10;
            }
        }
        return null;
    }
    _isLoopValid(loop) {
        if (!loop || loop.getEdges().length < 2) {
            return false;
        }
        const edges = loop.getEdges();
        for (let i = 0; i < edges.length; i++) {
            const edge1 = edges[i];
            for (let j = i + 1; j < edges.length; j++) {
                const edge2 = edges[j];
                if (this._doEdgesIntersect(edge1, edge2))
                    return false;
            }
        }
        const startVertex = edges[0].getStartVertex();
        const endVertex = edges[edges.length - 1].getEndVertex();
        if (startVertex !== endVertex)
            return false;
        return true;
    }
    _doEdgesIntersect(edge1, edge2) {
        if (!edge1 || !edge2) {
            return false;
        }
        const curve1 = edge1.getCurve();
        const curve2 = edge2.getCurve();
        if (!curve1 || !curve2) {
            const p1Start = edge1.getStartVertex().getPoint();
            const p1End = edge1.getEndVertex().getPoint();
            const p2Start = edge2.getStartVertex().getPoint();
            const p2End = edge2.getEndVertex().getPoint();
            if (this._doLineSegmentsIntersect(p1Start, p1End, p2Start, p2End)) {
                return true;
            } else {
                return false;
            }
        }
        const numSamples = 10;
        const start1 = curve1.getKnots()[0];
        const end1 = curve1.getKnots()[curve1.getKnots().length - 1];
        const step1 = (end1 - start1) / numSamples;
        const start2 = curve2.getKnots()[0];
        const end2 = curve2.getKnots()[curve2.getKnots().length - 1];
        const step2 = (end2 - start2) / numSamples;
        const points1 = [];
        const points2 = [];
        for (let i = 0; i <= numSamples; i++) {
            const u = start1 + i * step1;
            const point = curve1.evaluate(u);
            if (point) {
                points1.push(point);
            }
        }
        for (let i = 0; i <= numSamples; i++) {
            const u = start2 + i * step2;
            const point = curve2.evaluate(u);
            if (point) {
                points2.push(point);
            }
        }
        for (let i = 0; i < points1.length - 1; i++) {
            for (let j = 0; j < points2.length - 1; j++) {
                const p1Start = points1[i];
                const p1End = points1[i + 1];
                const p2Start = points2[j];
                const p2End = points2[j + 1];
                if (this._doLineSegmentsIntersect(p1Start, p1End, p2Start, p2End)) {
                    return true;
                }
            }
        }
        return false;
    }
    _calculateJacobianForIntersection(surface1, surface2, u1, v1, u2, v2) {
        const du1 = surface1.evaluatePartialDerivative(u1, v1, 1, 0);
        const dv1 = surface1.evaluatePartialDerivative(u1, v1, 0, 1);
        const du2 = surface2.evaluatePartialDerivative(u2, v2, 1, 0);
        const dv2 = surface2.evaluatePartialDerivative(u2, v2, 0, 1);
        if (!du1 || !dv1 || !du2 || !dv2) {
            return null;
        }
        return [
            [
                du1.getX(),
                du1.getY(),
                du1.getZ(),
                0,
                0,
                0
            ],
            [
                dv1.getX(),
                dv1.getY(),
                dv1.getZ(),
                0,
                0,
                0
            ],
            [
                0,
                0,
                0,
                -du2.getX(),
                -du2.getY(),
                -du2.getZ()
            ],
            [
                0,
                0,
                0,
                -dv2.getX(),
                -dv2.getY(),
                -dv2.getZ()
            ]
        ];
    }
    _transposeMultiply(matrixA, matrixB) {
        if (!matrixA || !matrixB || matrixA[0].length !== matrixB.length) {
            return null;
        }
        const rowsA = matrixA.length;
        const colsA = matrixA[0].length;
        const colsB = matrixB[0].length;
        const result = Array(colsA).fill(null).map(() => Array(colsB).fill(0));
        for (let i = 0; i < colsA; i++) {
            for (let j = 0; j < colsB; j++) {
                for (let k = 0; k < rowsA; k++) {
                    result[i][j] += matrixA[k][i] * matrixB[k][j];
                }
            }
        }
        return result;
    }
    _augmentMatrix(matrix, lambda) {
        if (!matrix) {
            return null;
        }
        const rows = matrix.length;
        const cols = matrix[0].length;
        if (rows !== cols) {
            return null;
        }
        const result = matrix.map(row => [...row]);
        for (let i = 0; i < rows; i++) {
            result[i][i] += lambda;
        }
        return result;
    }
    _solveLinearSystem(matrixA, matrixB) {
        if (!matrixA || !matrixB || matrixA.length !== matrixA[0].length || matrixA.length !== matrixB.length) {
            return null;
        }
        const n = matrixA.length;
        const augmentedMatrix = matrixA.map((row, i) => [
            ...row,
            ...matrixB[i]
        ]);
        for (let i = 0; i < n; i++) {
            let pivotRow = i;
            for (let j = i + 1; j < n; j++) {
                if (Math.abs(augmentedMatrix[j][i]) > Math.abs(augmentedMatrix[pivotRow][i])) {
                    pivotRow = j;
                }
            }
            if (pivotRow !== i) {
                [augmentedMatrix[i], augmentedMatrix[pivotRow]] = [
                    augmentedMatrix[pivotRow],
                    augmentedMatrix[i]
                ];
            }
            if (Math.abs(augmentedMatrix[i][i]) < 1e-8) {
                return null;
            }
            for (let j = 0; j < n; j++) {
                if (j !== i) {
                    const factor = augmentedMatrix[j][i] / augmentedMatrix[i][i];
                    for (let k = i; k < n + matrixB[0].length; k++) {
                        augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                    }
                }
            }
        }
        const result = matrixB.map((row, i) => []);
        for (let i = 0; i < n; i++) {
            for (let j = n; j < n + matrixB[0].length; j++) {
                result[i].push(augmentedMatrix[i][j] / augmentedMatrix[i][i]);
            }
        }
        return result;
    }
    _doLineSegmentsIntersect(p1Start, p1End, p2Start, p2End) {
        const x1 = p1Start.getX();
        const y1 = p1Start.getY();
        const x2 = p1End.getX();
        const y2 = p1End.getY();
        const x3 = p2Start.getX();
        const y3 = p2Start.getY();
        const x4 = p2End.getX();
        const y4 = p2End.getY();
        const det = (x2 - x1) * (y4 - y3) - (x4 - x3) * (y2 - y1);
        if (det === 0) {
            return false;
        }
        const t = ((x3 - x1) * (y4 - y3) - (x4 - x3) * (y3 - y1)) / det;
        const u = -((x1 - x2) * (y3 - y1) - (x3 - x1) * (y1 - y2)) / det;
        return t >= 0 && t <= 1 && u >= 0 && u <= 1;
    }
    _isSolidInside(solid1, solid2) {
        let allInside1 = true;
        for (const shell of solid1.getShells()) {
            for (const face of shell.getFaces()) {
                for (const loop of [
                        face.getOuterLoop(),
                        ...face.getInnerLoops()
                    ]) {
                    for (const edge of loop.getEdges()) {
                        const startPoint = edge.getStartVertex().getPoint();
                        const endPoint = edge.getEndVertex().getPoint();
                        const trimmer = new SurfaceTrimmer();
                        if (!trimmer.isPointInside(solid2.getFaces()[0].getSurface(), startPoint) || !trimmer.isPointInside(solid2.getFaces()[0].getSurface(), endPoint)) {
                            allInside1 = false;
                            break;
                        }
                    }
                    if (!allInside1)
                        break;
                }
                if (!allInside1)
                    break;
            }
            if (!allInside1)
                break;
        }
        if (allInside1)
            return true;
        let allInside2 = true;
        for (const shell of solid2.getShells()) {
            for (const face of shell.getFaces()) {
                for (const loop of [
                        face.getOuterLoop(),
                        ...face.getInnerLoops()
                    ]) {
                    for (const edge of loop.getEdges()) {
                        const startPoint = edge.getStartVertex().getPoint();
                        const endPoint = edge.getEndVertex().getPoint();
                        const trimmer = new SurfaceTrimmer();
                        if (!trimmer.isPointInside(solid1.getFaces()[0].getSurface(), startPoint) || !trimmer.isPointInside(solid1.getFaces()[0].getSurface(), endPoint)) {
                            allInside2 = false;
                            break;
                        }
                    }
                    if (!allInside2)
                        break;
                }
                if (!allInside2)
                    break;
            }
            if (!allInside2)
                break;
        }
        return allInside2;
    }
    // ... other methods ...
    _isCoplanar(surface1, surface2) {
        if (!surface1 || !surface2) {
            return false;
        }
        const uMin1 = surface1.getKnotsU()[0];
        const uMax1 = surface1.getKnotsU()[surface1.getKnotsU().length - 1];
        const vMin1 = surface1.getKnotsV()[0];
        const vMax1 = surface1.getKnotsV()[surface1.getKnotsV().length - 1];
        const uMin2 = surface2.getKnotsU()[0];
        const uMax2 = surface2.getKnotsU()[surface2.getKnotsU().length - 1];
        const vMin2 = surface2.getKnotsV()[0];
        const vMax2 = surface2.getKnotsV()[surface2.getKnotsV().length - 1];
        const numSamplesU = 5;
        const numSamplesV = 5;
        const deltaU1 = (uMax1 - uMin1) / numSamplesU;
        const deltaV1 = (vMax1 - vMin1) / numSamplesV;
        const tolerance = 0.001;
        let normal1 = surface1.normal((uMin1 + uMax1) / 2, (vMin1 + vMax1) / 2);
        let normal2 = surface2.normal((uMin2 + uMax2) / 2, (vMin2 + vMax2) / 2);
        if (!normal1 || !normal2)
            return false;
        const dotProduct = normal1.normalize().dot(normal2.normalize());
        if (Math.abs(Math.abs(dotProduct) - 1) > tolerance) {
            return false;
        }
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u1 = uMin1 + i * deltaU1;
                const v1 = vMin1 + j * deltaV1;
                const point1 = surface1.evaluate(u1, v1);
                if (!point1)
                    continue;
                const uv2 = this._findUVOnSurface(surface2, point1);
                if (!uv2)
                    return false;
            }
        }
        return true;
    }
    _handleCoplanarFaces(face1, face2, keepSide) {
        const surface1 = face1.getSurface();
        const surface2 = face2.getSurface();
        if (!surface1 || !surface2)
            return null;
        const uMin1 = surface1.getKnotsU()[0];
        const uMax1 = surface1.getKnotsU()[surface1.getKnotsU().length - 1];
        const vMin1 = surface1.getKnotsV()[0];
        const vMax1 = surface1.getKnotsV()[surface1.getKnotsV().length - 1];
        const uMin2 = surface2.getKnotsU()[0];
        const uMax2 = surface2.getKnotsU()[surface2.getKnotsU().length - 1];
        const vMin2 = surface2.getKnotsV()[0];
        const vMax2 = surface2.getKnotsV()[surface2.getKnotsV().length - 1];
        const numSamplesU = 10;
        const numSamplesV = 10;
        const deltaU1 = (uMax1 - uMin1) / numSamplesU;
        const deltaV1 = (vMax1 - vMin1) / numSamplesV;
        const points1 = [];
        const points2 = [];
        const face1OuterEdges = face1.getOuterLoop().getEdges();
        const face1InnerEdges = face1.getInnerLoops().flatMap(loop => loop.getEdges());
        const face2OuterEdges = face2.getOuterLoop().getEdges();
        const face2InnerEdges = face2.getInnerLoops().flatMap(loop => loop.getEdges());
        for (const edge of face1OuterEdges) {
            const startPoint = edge.getStartVertex().getPoint();
            const endPoint = edge.getEndVertex().getPoint();
            points1.push(startPoint, endPoint);
        }
        for (const edge of face1InnerEdges) {
            const startPoint = edge.getStartVertex().getPoint();
            const endPoint = edge.getEndVertex().getPoint();
            points1.push(startPoint, endPoint);
        }
        for (const edge of face2OuterEdges) {
            const startPoint = edge.getStartVertex().getPoint();
            const endPoint = edge.getEndVertex().getPoint();
            points2.push(startPoint, endPoint);
        }
        for (const edge of face2InnerEdges) {
            const startPoint = edge.getStartVertex().getPoint();
            const endPoint = edge.getEndVertex().getPoint();
            points2.push(startPoint, endPoint);
        }
        if (points1.length < 3 || points2.length < 3)
            return null;
        const face1Edges = this._createEdgesFromPoints(points1);
        const face2Edges = this._createEdgesFromPoints(points2);
        const face1Loop = new Loop(face1Edges);
        const face2Loop = new Loop(face2Edges);
        if (!this._isLoopValid(face1Loop) || !this._isLoopValid(face2Loop))
            return null;
        let trimmedSurface = null;
        let hasOverlap = false;
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u1 = uMin1 + i * deltaU1;
                const v1 = vMin1 + j * deltaV1;
                const point1 = surface1.evaluate(u1, v1);
                if (point1 && this._isPointInsideLoop(face1Loop, point1, surface1) && this._isPointInsideLoop(face2Loop, point1, surface2)) {
                    hasOverlap = true;
                }
            }
        }
        if (!hasOverlap)
            return null;
        trimmedSurface = this._trimCoplanarFaces(surface1, face1Loop, face2Loop, keepSide);
        if (trimmedSurface)
            return trimmedSurface;
        return null;
    }
    _createEdgesFromPoints(points) {
        if (!points || points.length < 2) {
            return null;
        }
        const edges = [];
        for (let i = 0; i < points.length - 1; i++) {
            const startPoint = points[i];
            const endPoint = points[i + 1];
            const startVertex = this._findOrCreateVertex(startPoint);
            const endVertex = this._findOrCreateVertex(endPoint);
            if (startVertex === endVertex) {
                continue;
            }
            const edgeCurve = new NurbsCurve(1, [
                startPoint,
                endPoint
            ], [
                0,
                1
            ], [
                1,
                1
            ]);
            const edge = new Edge(startVertex, endVertex, edgeCurve);
            edges.push(edge);
        }
        return edges;
    }
    _trimCoplanarFaces(surface, face1Loop, face2Loop, keepSide) {
        if (!surface || !face1Loop || !face2Loop)
            return null;
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        let newControlPoints = surface.getControlPoints().map(row => [...row]);
        let newWeights = surface.getWeights().map(row => [...row]);
        const numSamplesU = 20;
        const numSamplesV = 20;
        const uStep = (uMax - uMin) / numSamplesU;
        const vStep = (vMax - vMin) / numSamplesV;
        let hasValidControlPoints = false;
        for (let k = 0; k < newControlPoints.length; k++) {
            for (let l = 0; l < newControlPoints[k].length; l++) {
                const controlPoint = newControlPoints[k][l];
                const weight = newWeights[k][l];
                if (controlPoint === null || weight === null)
                    continue;
                const cpUv = this._findUVOnSurface(surface, controlPoint);
                if (cpUv) {
                    const isInsideFace1 = this._isPointInsideLoop(face1Loop, controlPoint, surface);
                    const isInsideFace2 = this._isPointInsideLoop(face2Loop, controlPoint, surface);
                    let isInside = false;
                    if (keepSide) {
                        isInside = isInsideFace1;
                    } else {
                        isInside = isInsideFace1 && !isInsideFace2;
                    }
                    if (!isInside) {
                        newControlPoints[k][l] = new Point3D(0, 0, 0);
                        newWeights[k][l] = 0;
                    } else {
                        hasValidControlPoints = true;
                    }
                }
            }
        }
        if (!hasValidControlPoints) {
            return null;
        }
        // Prevent creating a zero-area face
        return new NurbsSurface(surface.getDegreeU(), surface.getDegreeV(), newControlPoints, surface.getKnotsU(), surface.getKnotsV(), newWeights);
    }
    _isPointInsideLoop(loop, point, surface) {
        if (!loop || !point || !surface)
            return false;
        const uv = this._findUVOnSurface(surface, point);
        if (!uv)
            return false;
        let windingNumber = 0;
        const edges = loop.getEdges();
        for (let i = 0; i < edges.length; i++) {
            const edge = edges[i];
            const curve = edge.getCurve();
            if (!curve) {
                const p1 = edge.getStartVertex().getPoint();
                const p2 = edge.getEndVertex().getPoint();
                const uv1 = this._findUVOnSurface(surface, p1);
                const uv2 = this._findUVOnSurface(surface, p2);
                if (!uv1 || !uv2)
                    continue;
                if (uv1.v <= uv.v) {
                    if (uv2.v > uv.v) {
                        const vt = (uv.v - uv1.v) / (uv2.v - uv1.v);
                        if (uv.u < uv1.u + vt * (uv2.u - uv1.u)) {
                            windingNumber++;
                        }
                    }
                } else {
                    if (uv2.v <= uv.v) {
                        const vt = (uv.v - uv1.v) / (uv2.v - uv1.v);
                        if (uv.u < uv1.u + vt * (uv2.u - uv1.u)) {
                            windingNumber--;
                        }
                    }
                }
            } else {
                const numSamples = 10;
                const start = curve.getKnots()[0];
                const end = curve.getKnots()[curve.getKnots().length - 1];
                const step = (end - start) / numSamples;
                const points = [];
                for (let j = 0; j <= numSamples; j++) {
                    const u = start + j * step;
                    const p = curve.evaluate(u);
                    if (p)
                        points.push(p);
                }
                for (let j = 0; j < points.length - 1; j++) {
                    const p1 = points[j];
                    const p2 = points[j + 1];
                    const uv1 = this._findUVOnSurface(surface, p1);
                    const uv2 = this._findUVOnSurface(surface, p2);
                    if (!uv1 || !uv2)
                        continue;
                    if (uv1.v <= uv.v) {
                        if (uv2.v > uv.v) {
                            const vt = (uv.v - uv1.v) / (uv2.v - uv1.v);
                            if (uv.u < uv1.u + vt * (uv2.u - uv1.u)) {
                                windingNumber++;
                            }
                        }
                    } else {
                        if (uv2.v <= uv.v) {
                            const vt = (uv.v - uv1.v) / (uv2.v - uv1.v);
                            if (uv.u < uv1.u + vt * (uv2.u - uv1.u)) {
                                windingNumber--;
                            }
                        }
                    }
                }
            }
        }
        const isOuter = this._isOuterLoop(loop, surface);
        if (isOuter) {
            return windingNumber !== 0;
        } else {
            return windingNumber === 0;
        }
    }
    _getBoundingBox(surface, uMin, uMax, vMin, vMax) {
        const numSamplesU = 5;
        const numSamplesV = 5;
        const uStep = (uMax - uMin) / numSamplesU;
        const vStep = (vMax - vMin) / numSamplesV;
        let minX = Infinity, minY = Infinity, minZ = Infinity;
        let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const point = surface.evaluate(u, v);
                if (point) {
                    minX = Math.min(minX, point.getX());
                    minY = Math.min(minY, point.getY());
                    minZ = Math.min(minZ, point.getZ());
                    maxX = Math.max(maxX, point.getX());
                    maxY = Math.max(maxY, point.getY());
                    maxZ = Math.max(maxZ, point.getZ());
                }
            }
        }
        return {
            minX: minX,
            minY: minY,
            minZ: minZ,
            maxX: maxX,
            maxY: maxY,
            maxZ: maxZ
        };
    }
    _doBoundingBoxesOverlap(bbox1, bbox2) {
        return bbox1.maxX >= bbox2.minX && bbox1.minX <= bbox2.maxX && bbox1.maxY >= bbox2.minY && bbox1.minY <= bbox2.maxY && bbox1.maxZ >= bbox2.minZ && bbox1.minZ <= bbox2.maxZ;
    }
    _trimSurfaceWithUVs(surface, uvs, keepSide) {
        if (!surface || !uvs || uvs.length < 2) {
            return surface;
        }
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        let newControlPoints = surface.getControlPoints().map(row => [...row]);
        let newWeights = surface.getWeights().map(row => [...row]);
        const numSamplesU = 20;
        const numSamplesV = 20;
        const uStep = (uMax - uMin) / numSamplesU;
        const vStep = (vMax - vMin) / numSamplesV;
        let hasValidControlPoints = false;
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const point = surface.evaluate(u, v);
                if (!point)
                    continue;
                const isInside = this._isPointInsideTrim(uvs, {
                    u: u,
                    v: v
                }, keepSide);
                if (!isInside) {
                    for (let k = 0; k < newControlPoints.length; k++) {
                        for (let l = 0; l < newControlPoints[k].length; l++) {
                            const controlPoint = newControlPoints[k][l];
                            const weight = newWeights[k][l];
                            if (controlPoint === null || weight === null)
                                continue;
                            const cpUv = this._findUVOnSurface(surface, controlPoint);
                            if (cpUv) {
                                const dist = Math.sqrt((cpUv.u - u) * (cpUv.u - u) + (cpUv.v - v) * (cpUv.v - v));
                                if (dist < 0.1) {
                                    newControlPoints[k][l] = new Point3D(0, 0, 0);
                                    newWeights[k][l] = 0;
                                }
                            }
                        }
                    }
                } else {
                    hasValidControlPoints = true;
                }
            }
        }
        if (!hasValidControlPoints) {
            return null;
        }
        // Prevent creating a zero-area face
        return new NurbsSurface(surface.getDegreeU(), surface.getDegreeV(), newControlPoints, surface.getKnotsU(), surface.getKnotsV(), newWeights);
    }
    _isOuterLoop(loop, surface) {
        if (!loop || !surface)
            return false;
        const edges = loop.getEdges();
        if (edges.length === 0)
            return false;
        let area = 0;
        for (const edge of edges) {
            const curve = edge.getCurve();
            if (curve) {
                const numSamples = 10;
                const start = curve.getKnots()[0];
                const end = curve.getKnots()[curve.getKnots().length - 1];
                const step = (end - start) / numSamples;
                const points = [];
                for (let j = 0; j <= numSamples; j++) {
                    const u = start + j * step;
                    const p = curve.evaluate(u);
                    if (p)
                        points.push(p);
                }
                for (let i = 0; i < points.length - 1; i++) {
                    const p1 = points[i];
                    const p2 = points[i + 1];
                    area += (p1.getX() * p2.getY() - p2.getX() * p1.getY()) * 0.5;
                }
            } else {
                const p1 = edge.getStartVertex().getPoint();
                const p2 = edge.getEndVertex().getPoint();
                area += (p1.getX() * p2.getY() - p2.getX() * p1.getY()) * 0.5;
            }
        }
        if (area >= 0) {
            return true;
        }
        return false;
    }
    _isFaceLargeEnough(face) {
        if (!face || !face.getSurface()) {
            return false;
        }
        const surface = face.getSurface();
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        const numSamplesU = 10;
        const numSamplesV = 10;
        const deltaU = (uMax - uMin) / numSamplesU;
        const deltaV = (vMax - vMin) / numSamplesV;
        let area = 0;
        for (let i = 0; i < numSamplesU; i++) {
            for (let j = 0; j < numSamplesV; j++) {
                const u = uMin + i * deltaU;
                const v = vMin + j * deltaV;
                const nextU = u + deltaU;
                const nextV = v + deltaV;
                const p0 = surface.evaluate(u, v);
                const p1 = surface.evaluate(nextU, v);
                const p2 = surface.evaluate(u, nextV);
                if (p0 && p1 && p2) {
                    const v1 = new Vector3D(p1.getX() - p0.getX(), p1.getY() - p0.getY(), p1.getZ() - p0.getZ());
                    const v2 = new Vector3D(p2.getX() - p0.getX(), p2.getY() - p0.getY(), p2.getZ() - p0.getZ());
                    area += v1.cross(v2).magnitude() * 0.5;
                }
            }
        }
        const bbox = this._getBoundingBox(surface, uMin, uMax, vMin, vMax);
        const width = bbox.maxX - bbox.minX;
        const height = bbox.maxY - bbox.minY;
        const depth = bbox.maxZ - bbox.minZ;
        if (area > 0.0001 && (width > 0.0001 || height > 0.0001 || depth > 0.0001))
            return true;
        return false;
    }
    _isFaceInsideLoop(loop1, loop2, surface) {
        if (!loop1 || !loop2 || !surface) {
            return false;
        }
        const trimmer = new SurfaceTrimmer();
        for (const edge of loop2.getEdges()) {
            const startPoint = edge.getStartVertex().getPoint();
            const endPoint = edge.getEndVertex().getPoint();
            if (!trimmer.isPointInside(surface, startPoint) || !trimmer.isPointInside(surface, endPoint)) {
                return false;
            }
        }
        return true;
    }
    _isCurveWithinBounds(curve, surface) {
        if (!curve || !surface) {
            return false;
        }
        const numSamples = 10;
        const start = curve.getKnots()[0];
        const end = curve.getKnots()[curve.getKnots().length - 1];
        const step = (end - start) / numSamples;
        for (let i = 0; i <= numSamples; i++) {
            const u = start + i * step;
            const point = curve.evaluate(u);
            if (!point)
                continue;
            const uv = this._findUVOnSurface(surface, point);
            if (!uv)
                return false;
            const uMin = surface.getKnotsU()[0];
            const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
            const vMin = surface.getKnotsV()[0];
            const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
            const tolerance = 0.001;
            if (uv.u < uMin - tolerance || uv.u > uMax + tolerance || uv.v < vMin - tolerance || uv.v > vMax + tolerance) {
                return false;
            }
        }
        return true;
    }
    _splitSelfIntersectingSurface(surface) {
        if (!surface) {
            return null;
        }
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        const numSamplesU = 10;
        const numSamplesV = 10;
        const uStep = (uMax - uMin) / numSamplesU;
        const vStep = (vMax - vMin) / numSamplesV;
        const points = [];
        for (let i = 0; i <= numSamplesU; i++) {
            for (let j = 0; j <= numSamplesV; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const point = surface.evaluate(u, v);
                if (point)
                    points.push({
                        u,
                        v,
                        point
                    });
            }
        }
        const intersections = [];
        for (let i = 0; i < points.length; i++) {
            for (let j = i + 1; j < points.length; j++) {
                const point1 = points[i].point;
                const point2 = points[j].point;
                const dist = point1.distanceTo(point2);
                if (dist < 0.001) {
                    const u1 = points[i].u;
                    const v1 = points[i].v;
                    const u2 = points[j].u;
                    const v2 = points[j].v;
                    if (Math.abs(u1 - u2) > 0.001 || Math.abs(v1 - v2) > 0.001) {
                        intersections.push({
                            u: u1,
                            v: v1,
                            u2: u2,
                            v2: v2
                        });
                    }
                }
            }
        }
        if (intersections.length === 0) {
            return [surface];
        }
        let splitSurfaces = [surface];
        const trimmer = new SurfaceTrimmer();
        for (const intersection of intersections) {
            let newSplitSurfaces = [];
            for (const surface of splitSurfaces) {
                const refinedIntersection = this._refineIntersectionPoint(surface, surface, {
                    u1: intersection.u,
                    v1: intersection.v,
                    u2: intersection.u2,
                    v2: intersection.v2
                }, 0.001, 20);
                if (refinedIntersection) {
                    const splitResult = trimmer.splitSurface(surface, new Intersection(new NurbsCurve(1, [
                        surface.evaluate(refinedIntersection.u1, refinedIntersection.v1),
                        surface.evaluate(refinedIntersection.u2, refinedIntersection.v2)
                    ], [
                        0,
                        1
                    ], [
                        1,
                        1
                    ]), surface, surface));
                    if (splitResult && splitResult.length > 0) {
                        newSplitSurfaces.push(...splitResult);
                    } else {
                        newSplitSurfaces.push(surface);
                    }
                } else {
                    newSplitSurfaces.push(surface);
                }
            }
            splitSurfaces = newSplitSurfaces;
        }
        return splitSurfaces;
    }
    isPointInside(surface, point) {
        if (!surface || !point) {
            return false;
        }
        const uv = this._findUVOnSurface(surface, point);
        if (!uv) {
            return false;
        }
        const uvs = [];
        const numSamples = 50;
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        const uStep = (uMax - uMin) / numSamples;
        const vStep = (vMax - vMin) / numSamples;
        for (let i = 0; i <= numSamples; i++) {
            for (let j = 0; j <= numSamples; j++) {
                const u = uMin + i * uStep;
                const v = vMin + j * vStep;
                const p = surface.evaluate(u, v);
                if (p) {
                    const dist = p.distanceTo(point);
                    if (dist < 0.001) {
                        uvs.push({
                            u: u,
                            v: v
                        });
                    }
                }
            }
        }
        if (uvs.length === 0) {
            return false;
        }
        const testUv = uvs[0];
        if (!testUv)
            return false;
        return this._isPointInsideTrim(uvs, testUv, true);
    }
    _findUVOnSurface(surface, point) {
        const uMin = surface.getKnotsU()[0];
        const uMax = surface.getKnotsU()[surface.getKnotsU().length - 1];
        const vMin = surface.getKnotsV()[0];
        const vMax = surface.getKnotsV()[surface.getKnotsV().length - 1];
        let u = (uMin + uMax) / 2;
        let v = (vMin + vMax) / 2;
        const tolerance = 0.001;
        const maxIterations = 20;
        for (let i = 0; i < maxIterations; i++) {
            const surfacePoint = surface.evaluate(u, v);
            if (!surfacePoint) {
                return null;
            }
            const diff = new Vector3D(surfacePoint.getX() - point.getX(), surfacePoint.getY() - point.getY(), surfacePoint.getZ() - point.getZ());
            if (diff.magnitude() < tolerance) {
                if (u >= uMin && u <= uMax && v >= vMin && v <= vMax)
                    return {
                        u: u,
                        v: v
                    };
                else {
                    return null;
                }
            }
            const jacobian = this._calculateJacobian(surface, u, v);
            if (!jacobian) {
                return null;
            }
            const invJacobian = this._inverse2x2(jacobian);
            if (!invJacobian) {
                return null;
            }
            const deltaU = invJacobian[0][0] * diff.getX() + invJacobian[0][1] * diff.getY() + invJacobian[0][2] * diff.getZ();
            const deltaV = invJacobian[1][0] * diff.getX() + invJacobian[1][1] * diff.getY() + invJacobian[1][2] * diff.getZ();
            u -= deltaU;
            v -= deltaV;
            u = Math.max(uMin, Math.min(u, uMax));
            v = Math.max(vMin, Math.min(v, vMax));
        }
        return null;
    }
    _inverse2x2(matrix) {
        const det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        if (Math.abs(det) < 1e-8) {
            return null;
        }
        const invDet = 1 / det;
        return [
            [
                matrix[1][1] * invDet,
                -matrix[0][1] * invDet,
                0
            ],
            [
                -matrix[1][0] * invDet,
                matrix[0][0] * invDet,
                0
            ]
        ];
    }
}