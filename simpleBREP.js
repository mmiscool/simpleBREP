
export class Point {
    constructor(x, y, z) {
        this.x = x || 0;
        // Default to 0 if not provided
        this.y = y || 0;
        this.z = z || 0;
    }
    add(point) {
        return new Point(this.x + point.x, this.y + point.y, this.z + point.z);
    }
    subtract(point) {
        return new Point(this.x - point.x, this.y - point.y, this.z - point.z);
    }
    scale(scalar) {
        return new Point(this.x * scalar, this.y * scalar, this.z * scalar);
    }
    distanceTo(point) {
        const dx = this.x - point.x;
        const dy = this.y - point.y;
        const dz = this.z - point.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
    equals(point, tolerance = 0.000001) {
        return this.distanceTo(point) < tolerance;
    }
    clone() {
        return new Point(this.x, this.y, this.z);
    }
    toString() {
        return `Point(${ this.x }, ${ this.y }, ${ this.z })`;
    }
}
export class Vector {
    constructor(x, y, z) {
        this.x = x || 0;
        // Default to 0 if not provided
        this.y = y || 0;
        this.z = z || 0;
    }
    add(vector) {
        return new Vector(this.x + vector.x, this.y + vector.y, this.z + vector.z);
    }
    subtract(vector) {
        return new Vector(this.x - vector.x, this.y - vector.y, this.z - vector.z);
    }
    scale(scalar) {
        return new Vector(this.x * scalar, this.y * scalar, this.z * scalar);
    }
    magnitude() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }
    normalize() {
        const mag = this.magnitude();
        if (mag === 0) {
            throw new Error('Cannot normalize a zero-length vector');
        }
        return this.scale(1 / mag);
    }
    dot(vector) {
        return this.x * vector.x + this.y * vector.y + this.z * vector.z;
    }
    cross(vector) {
        return new Vector(this.y * vector.z - this.z * vector.y, this.z * vector.x - this.x * vector.z, this.x * vector.y - this.y * vector.x);
    }
    equals(vector, tolerance = 0.000001) {
        const diff = this.subtract(vector);
        return diff.magnitude() < tolerance;
    }
    clone() {
        return new Vector(this.x, this.y, this.z);
    }
    toString() {
        return `Vector(${ this.x }, ${ this.y }, ${ this.z })`;
    }
}
export class NURBSCurve {
    constructor(degree, controlPoints, knots, weights) {
        if (!Number.isInteger(degree) || degree < 1) {
            throw new Error('Degree must be a positive integer');
        }
        if (!Array.isArray(controlPoints) || controlPoints.length < degree + 1) {
            throw new Error('Must have at least degree + 1 control points');
        }
        if (!Array.isArray(knots) || knots.length !== controlPoints.length + degree + 1) {
            throw new Error('Knot vector length must be controlPoints.length + degree + 1');
        }
        if (!Array.isArray(weights) || weights.length !== controlPoints.length) {
            throw new Error('Weights array must match controlPoints length');
        }
        if (!knots.every((k, i) => i === 0 || k >= knots[i - 1])) {
            throw new Error('Knots must be non-decreasing');
        }
        this.degree = degree;
        this.controlPoints = controlPoints.map(p => p.clone());
        this.knots = [...knots];
        this.weights = [...weights];
        this.n = controlPoints.length - 1;
        this._basisCache = new Map();
        this._derivativeCache = new Map();
    }
    evaluate(t) {
        const span = this._findKnotSpan(t);
        const basis = this._getCachedBasis(span, t);
        let x = 0, y = 0, z = 0, w = 0;
        for (let i = 0; i <= this.degree; i++) {
            const idx = span - this.degree + i;
            const b = basis[i] * this.weights[idx];
            x += this.controlPoints[idx].x * b;
            y += this.controlPoints[idx].y * b;
            z += this.controlPoints[idx].z * b;
            w += b;
        }
        return new Point(x / w, y / w, z / w);
    }
    _findKnotSpan(t) {
        if (t < this.knots[0] || t > this.knots[this.knots.length - 1]) {
            throw new Error(`Parameter t=${ t } is outside knot range [${ this.knots[0] }, ${ this.knots[this.knots.length - 1] }]`);
        }
        if (t === this.knots[this.knots.length - 1]) {
            return this.n;
        }
        let low = this.degree;
        let high = this.n + 1;
        let mid = Math.floor((low + high) / 2);
        while (t < this.knots[mid] || t >= this.knots[mid + 1]) {
            if (t < this.knots[mid])
                high = mid;
            else
                low = mid;
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    _basisFunctions(span, t) {
        const N = [1];
        const left = [];
        const right = [];
        for (let j = 1; j <= this.degree; j++) {
            left[j] = t - this.knots[span + 1 - j];
            right[j] = this.knots[span + j] - t;
            let saved = 0;
            for (let r = 0; r < j; r++) {
                const temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
        return N;
    }
    getDomain() {
        return [
            this.knots[this.degree],
            this.knots[this.n + 1]
        ];
    }
    toString() {
        return `NURBSCurve(degree=${ this.degree }, n=${ this.n }, knots=[${ this.knots.join(', ') }])`;
    }
    derivative(t, order = 1) {
        if (order < 1 || !Number.isInteger(order)) {
            throw new Error('Derivative order must be a positive integer');
        }
        if (order > this.degree) {
            return new Vector(0, 0, 0);
        }
        const span = this._findKnotSpan(t);
        const ders = this._getCachedDerivatives(span, t, order);
        const Aders = [];
        const wders = [];
        for (let k = 0; k <= order; k++) {
            let x = 0, y = 0, z = 0, w = 0;
            for (let i = 0; i <= this.degree; i++) {
                const idx = span - this.degree + i;
                const b = ders[k][i] * this.weights[idx];
                x += this.controlPoints[idx].x * b;
                y += this.controlPoints[idx].y * b;
                z += this.controlPoints[idx].z * b;
                w += ders[k][i] * this.weights[idx];
            }
            Aders[k] = new Point(x, y, z);
            wders[k] = w;
        }
        const result = new Vector(0, 0, 0);
        const w = wders[0];
        const wInv = 1 / w;
        const wInv2 = wInv * wInv;
        if (order === 1) {
            const A1 = Aders[1];
            const w1 = wders[1];
            result.x = (A1.x - w1 * Aders[0].x) * wInv2;
            result.y = (A1.y - w1 * Aders[0].y) * wInv2;
            result.z = (A1.z - w1 * Aders[0].z) * wInv2;
        } else {
            for (let k = 0; k <= order; k++) {
                const coef = this._binomial(order, k);
                const sign = (order - k) % 2 === 0 ? 1 : -1;
                let term = Aders[k].clone();
                for (let i = 0; i < order - k; i++) {
                    term = term.scale(wInv);
                }
                const wTerm = wders[order - k] * coef * sign;
                result.x += (term.x - Aders[0].x * wTerm) * wInv2;
                result.y += (term.y - Aders[0].y * wTerm) * wInv2;
                result.z += (term.z - Aders[0].z * wTerm) * wInv2;
            }
        }
        return result;
    }
    _basisDerivatives(span, t, n) {
        const ders = Array(n + 1).fill().map(() => Array(this.degree + 1).fill(0));
        const ndu = Array(this.degree + 1).fill().map(() => Array(this.degree + 1).fill(0));
        const left = Array(this.degree + 1).fill(0);
        const right = Array(this.degree + 1).fill(0);
        ndu[0][0] = 1;
        for (let j = 1; j <= this.degree; j++) {
            left[j] = t - this.knots[span + 1 - j];
            right[j] = this.knots[span + j] - t;
            let saved = 0;
            for (let r = 0; r < j; r++) {
                ndu[j][r] = right[r + 1] + left[j - r];
                const temp = ndu[r][j - 1] / ndu[j][r];
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }
        for (let j = 0; j <= this.degree; j++) {
            ders[0][j] = ndu[j][this.degree];
        }
        for (let r = 0; r <= this.degree; r++) {
            let s1 = 0, s2 = 1;
            const a = Array(this.degree + 1).fill().map(() => Array(this.degree + 1).fill(0));
            a[0][0] = 1;
            for (let k = 1; k <= n; k++) {
                let d = 0;
                const rk = r - k;
                const pk = this.degree - k;
                if (r >= k) {
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                    d = a[s2][0] * ndu[rk][pk];
                }
                const j1 = rk >= -1 ? 1 : -rk;
                const j2 = r - 1 <= pk ? k - 1 : this.degree - r;
                for (let j = j1; j <= j2; j++) {
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                    d += a[s2][j] * ndu[rk + j][pk];
                }
                if (r <= pk) {
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                    d += a[s2][k] * ndu[r][pk];
                }
                ders[k][r] = d;
                [s1, s2] = [
                    s2,
                    s1
                ];
            }
        }
        let r = this.degree;
        for (let k = 1; k <= n; k++) {
            for (let j = 0; j <= this.degree; j++) {
                ders[k][j] *= r;
            }
            r *= this.degree - k;
        }
        return ders;
    }
    closestParameter(point, tolerance = 0.000001, maxIterations = 100) {
        const [tMin, tMax] = this.getDomain();
        let t = (tMin + tMax) / 2;
        let iteration = 0;
        while (iteration < maxIterations) {
            const P = this.evaluate(t);
            const dP = this.derivative(t, 1);
            const V = P.subtract(point);
            const dot = V.dot(dP);
            const dPmag2 = dP.dot(dP);
            if (dPmag2 < tolerance)
                break;
            const deltaT = -dot / dPmag2;
            t += deltaT;
            if (t < tMin)
                t = tMin;
            if (t > tMax)
                t = tMax;
            if (Math.abs(deltaT) < tolerance)
                break;
            iteration++;
        }
        return t;
    }
    length(samples = 100) {
        const [tMin, tMax] = this.getDomain();
        const dt = (tMax - tMin) / (samples - 1);
        let totalLength = 0;
        let prevPoint = this.evaluate(tMin);
        for (let i = 1; i < samples; i++) {
            const t = tMin + i * dt;
            const currPoint = this.evaluate(t);
            totalLength += prevPoint.distanceTo(currPoint);
            prevPoint = currPoint;
        }
        return totalLength;
    }
    _binomial(n, k) {
        if (k < 0 || k > n)
            return 0;
        let result = 1;
        for (let i = 0; i < k; i++) {
            result *= (n - i) / (i + 1);
        }
        return result;
    }
    split(t) {
        const [tMin, tMax] = this.getDomain();
        if (t <= tMin || t >= tMax) {
            throw new Error(`Split parameter t=${ t } must be within domain (${ tMin }, ${ tMax })`);
        }
        let span = this._findKnotSpan(t);
        let multiplicity = 0;
        for (let i = span; i >= 0 && Math.abs(this.knots[i] - t) < 0.000001; i--) {
            multiplicity++;
        }
        const insertsNeeded = this.degree - multiplicity;
        let newKnots = [...this.knots];
        let newPoints = this.controlPoints.map(p => p.clone());
        let newWeights = [...this.weights];
        for (let i = 0; i < insertsNeeded; i++) {
            ({
                knots: newKnots,
                points: newPoints,
                weights: newWeights
            } = this._insertKnot(t, newKnots, newPoints, newWeights));
        }
        let splitIdx = 0;
        for (let i = 0; i < newKnots.length; i++) {
            if (newKnots[i] >= t - 0.000001) {
                splitIdx = i - this.degree;
                break;
            }
        }
        const knots1 = newKnots.slice(0, splitIdx + this.degree + 1);
        const points1 = newPoints.slice(0, splitIdx + 1);
        const weights1 = newWeights.slice(0, splitIdx + 1);
        const curve1 = new NURBSCurve(this.degree, points1, knots1, weights1);
        const knots2 = newKnots.slice(splitIdx);
        const points2 = newPoints.slice(splitIdx);
        const weights2 = newWeights.slice(splitIdx);
        const curve2 = new NURBSCurve(this.degree, points2, knots2, weights2);
        return [
            curve1,
            curve2
        ];
    }
    _insertKnot(u, knots, points, weights) {
        const span = this._findKnotSpan(u);
        const newKnots = [
            ...knots.slice(0, span + 1),
            u,
            ...knots.slice(span + 1)
        ];
        const newPoints = [];
        const newWeights = [];
        for (let i = 0; i <= span - this.degree; i++) {
            newPoints[i] = points[i].clone();
            newWeights[i] = weights[i];
        }
        for (let i = span - this.degree + 1; i < points.length; i++) {
            newPoints[i + 1] = points[i].clone();
            newWeights[i + 1] = weights[i];
        }
        const Q = [];
        const W = [];
        for (let i = 0; i <= this.degree; i++) {
            const idx = span - this.degree + i;
            const alpha = (u - knots[idx]) / (knots[idx + this.degree] - knots[idx]);
            Q[i] = points[idx - 1].scale(1 - alpha).add(points[idx].scale(alpha));
            W[i] = weights[idx - 1] * (1 - alpha) + weights[idx] * alpha;
        }
        for (let i = 0; i <= this.degree; i++) {
            newPoints[span - this.degree + i] = Q[i];
            newWeights[span - this.degree + i] = W[i];
        }
        return {
            knots: newKnots,
            points: newPoints,
            weights: newWeights
        };
    }
    _getCachedBasis(span, t) {
        const key = `${ span },${ t.toFixed(6) }`;
        if (!this._basisCache.has(key)) {
            this._basisCache.set(key, this._basisFunctions(span, t));
        }
        return this._basisCache.get(key);
    }
    _getCachedDerivatives(span, t, order) {
        const key = `${ span },${ t.toFixed(6) },${ order }`;
        if (!this._derivativeCache.has(key)) {
            this._derivativeCache.set(key, this._basisDerivatives(span, t, order));
        }
        return this._derivativeCache.get(key);
    }
    trim(t1, t2) {
        const [tMin, tMax] = this.getDomain();
        if (t1 < tMin || t2 > tMax || t1 >= t2) {
            throw new Error(`Trim parameters t1=${ t1 }, t2=${ t2 } must satisfy ${ tMin } <= t1 < t2 <= ${ tMax }`);
        }
        const [curveA, curveB] = this.split(t1);
        const [curveC, curveD] = curveB.split(t2 - t1);
        return curveC;
    }
    intersect(otherCurve, tolerance = 0.000001, maxIterations = 100) {
        const [tMin1, tMax1] = this.getDomain();
        const [tMin2, tMax2] = otherCurve.getDomain();
        const intersections = [];
        const subdivide = (t1a, t1b, t2a, t2b, depth = 0) => {
            if (depth > 10)
                return;
            const aabb1 = this._computeAABB(t1a, t1b);
            const aabb2 = otherCurve._computeAABB(t2a, t2b);
            if (!this._aabbIntersect(aabb1, aabb2))
                return;
            const p1a = this.evaluate(t1a);
            const p1b = this.evaluate(t1b);
            const p2a = otherCurve.evaluate(t2a);
            const p2b = otherCurve.evaluate(t2b);
            const d1 = p1a.distanceTo(p1b);
            const d2 = p2a.distanceTo(p2b);
            if (d1 < tolerance && d2 < tolerance) {
                const t1 = (t1a + t1b) / 2;
                const t2 = (t2a + t2b) / 2;
                if (this.evaluate(t1).distanceTo(otherCurve.evaluate(t2)) < tolerance) {
                    intersections.push({
                        t1,
                        t2
                    });
                }
                return;
            }
            const t1m = (t1a + t1b) / 2;
            const t2m = (t2a + t2b) / 2;
            subdivide(t1a, t1m, t2a, t2m, depth + 1);
            subdivide(t1a, t1m, t2m, t2b, depth + 1);
            subdivide(t1m, t1b, t2a, t2m, depth + 1);
            subdivide(t1m, t1b, t2m, t2b, depth + 1);
        };
        subdivide(tMin1, tMax1, tMin2, tMax2);
        const refined = [];
        for (let {t1, t2} of intersections) {
            let t1Refined = t1, t2Refined = t2;
            for (let i = 0; i < maxIterations; i++) {
                const P1 = this.evaluate(t1Refined);
                const P2 = otherCurve.evaluate(t2Refined);
                const V = P1.subtract(P2);
                const dC1 = this.derivative(t1Refined, 1);
                const dC2 = otherCurve.derivative(t2Refined, 1);
                const F = [
                    V.x,
                    V.y,
                    V.z
                ];
                const J = [
                    [
                        dC1.x,
                        -dC2.x
                    ],
                    [
                        dC1.y,
                        -dC2.y
                    ],
                    [
                        dC1.z,
                        -dC2.z
                    ]
                ];
                const JTJ = [
                    [
                        J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0],
                        J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]
                    ],
                    [
                        J[0][1] * J[0][0] + J[1][1] * J[1][0] + J[2][1] * J[2][0],
                        J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]
                    ]
                ];
                const JTF = [
                    J[0][0] * F[0] + J[1][0] * F[1] + J[2][0] * F[2],
                    J[0][1] * F[0] + J[1][1] * F[1] + J[2][1] * F[2]
                ];
                const det = JTJ[0][0] * JTJ[1][1] - JTJ[0][1] * JTJ[1][0];
                if (Math.abs(det) < tolerance)
                    break;
                const deltaT1 = (JTJ[1][1] * JTF[0] - JTJ[0][1] * JTF[1]) / det;
                const deltaT2 = (-JTJ[1][0] * JTF[0] + JTJ[0][0] * JTF[1]) / det;
                t1Refined -= deltaT1;
                t2Refined -= deltaT2;
                t1Refined = Math.max(tMin1, Math.min(tMax1, t1Refined));
                t2Refined = Math.max(tMin2, Math.min(tMax2, t2Refined));
                if (Math.abs(deltaT1) < tolerance && Math.abs(deltaT2) < tolerance) {
                    if (this.evaluate(t1Refined).distanceTo(otherCurve.evaluate(t2Refined)) < tolerance) {
                        refined.push({
                            t1: t1Refined,
                            t2: t2Refined
                        });
                    }
                    break;
                }
            }
        }
        return refined;
    }
    clearCache() {
        this._basisCache.clear();
        this._derivativeCache.clear();
    }
    _computeAABB(tStart, tEnd, samples = 5) {
        const dt = (tEnd - tStart) / (samples - 1);
        let minX = Infinity, minY = Infinity, minZ = Infinity;
        let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
        for (let i = 0; i < samples; i++) {
            const t = tStart + i * dt;
            const p = this.evaluate(t);
            minX = Math.min(minX, p.x);
            minY = Math.min(minY, p.y);
            minZ = Math.min(minZ, p.z);
            maxX = Math.max(maxX, p.x);
            maxY = Math.max(maxY, p.y);
            maxZ = Math.max(maxZ, p.z);
        }
        return {
            min: new Point(minX, minY, minZ),
            max: new Point(maxX, maxY, maxZ)
        };
    }
    _aabbIntersect(aabb1, aabb2) {
        return aabb1.min.x <= aabb2.max.x && aabb1.max.x >= aabb2.min.x && aabb1.min.y <= aabb2.max.y && aabb1.max.y >= aabb2.min.y && aabb1.min.z <= aabb2.max.z && aabb1.max.z >= aabb2.min.z;
    }
    elevateDegree() {
        const newDegree = this.degree + 1;
        const n = this.n;
        const newControlPoints = Array(n + 2);
        const newWeights = Array(n + 2);
        const newKnots = Array(n + newDegree + 2);
        // Copy boundary knots
        for (let i = 0; i <= this.degree; i++) {
            newKnots[i] = this.knots[0];
            newKnots[newKnots.length - 1 - i] = this.knots[this.knots.length - 1];
        }
        // Compute new control points and weights
        newControlPoints[0] = this.controlPoints[0].clone();
        newWeights[0] = this.weights[0];
        for (let i = 1; i <= n + 1; i++) {
            const alpha = i / (n + 1);
            newControlPoints[i] = this.controlPoints[i - 1].scale(1 - alpha).add(this.controlPoints[Math.min(i, n)].scale(alpha));
            newWeights[i] = this.weights[i - 1] * (1 - alpha) + this.weights[Math.min(i, n)] * alpha;
        }
        // Insert interior knots (uniform distribution for simplicity)
        for (let i = this.degree + 1; i <= n + 1; i++) {
            newKnots[i] = (i - this.degree) / (n - this.degree + 2);
        }
        return new NURBSCurve(newDegree, newControlPoints, newKnots, newWeights);
    }
    static fitPoints(points, degree = 3) {
        if (!Array.isArray(points) || points.length < degree + 1) {
            throw new Error('At least degree + 1 points are required for fitting');
        }
        const n = points.length - 1;
        const knots = Array(n + degree + 2).fill(0);
        const weights = Array(points.length).fill(1);
        // Generate uniform knots
        for (let i = degree + 1; i <= n; i++) {
            knots[i] = (i - degree) / (n - degree + 1);
        }
        for (let i = n + 1; i <= n + degree + 1; i++) {
            knots[i] = 1;
        }
        // Simple interpolation (centripetal parameterization)
        const params = [];
        let totalLength = 0;
        for (let i = 1; i < points.length; i++) {
            totalLength += Math.sqrt(points[i].distanceTo(points[i - 1]));
        }
        let cumulative = 0;
        params.push(0);
        for (let i = 1; i < points.length - 1; i++) {
            cumulative += Math.sqrt(points[i].distanceTo(points[i - 1]));
            params.push(cumulative / totalLength);
        }
        params.push(1);
        // For simplicity, use input points as control points (refinement can be added later)
        const controlPoints = points.map(p => p.clone());
        return new NURBSCurve(degree, controlPoints, knots, weights);
    }
    curvature(t) {
        const d1 = this.derivative(t, 1);
        const d2 = this.derivative(t, 2);
        const d1Mag = d1.magnitude();
        if (d1Mag < 0.000001) {
            return 0;
        }
        // Avoid division by zero
        const numerator = d1.cross(d2).magnitude();
        const denominator = d1Mag * d1Mag * d1Mag;
        return numerator / denominator;
    }
    reduceDegree(tolerance = 0.000001) {
        if (this.degree <= 1) {
            throw new Error('Cannot reduce degree below 1');
        }
        const newDegree = this.degree - 1;
        const n = this.n;
        const newControlPoints = Array(n).fill(null);
        const newWeights = Array(n).fill(0);
        const newKnots = Array(n + newDegree + 2).fill(0);
        // Copy boundary knots
        for (let i = 0; i <= newDegree; i++) {
            newKnots[i] = this.knots[0];
            newKnots[n + newDegree + 1 - i] = this.knots[this.knots.length - 1];
        }
        for (let i = newDegree + 1; i <= n; i++) {
            newKnots[i] = (i - newDegree) / (n - newDegree + 1);
        }
        // Approximate control points by least-squares fitting
        const samples = Math.max(n * 2, 10);
        const domain = this.getDomain();
        const dt = (domain[1] - domain[0]) / (samples - 1);
        const samplePoints = [];
        for (let i = 0; i < samples; i++) {
            const t = domain[0] + i * dt;
            samplePoints.push(this.evaluate(t));
        }
        // Fit a lower-degree curve to the sampled points
        const tempCurve = NURBSCurve.fitPoints(samplePoints, newDegree);
        for (let i = 0; i < n; i++) {
            newControlPoints[i] = tempCurve.controlPoints[i] || tempCurve.controlPoints[tempCurve.controlPoints.length - 1];
            newWeights[i] = tempCurve.weights[i] || tempCurve.weights[tempCurve.weights.length - 1];
        }
        const reducedCurve = new NURBSCurve(newDegree, newControlPoints, newKnots, newWeights);
        // Verify approximation within tolerance
        for (let i = 0; i < samples; i++) {
            const t = domain[0] + i * dt;
            const originalPoint = this.evaluate(t);
            const reducedPoint = reducedCurve.evaluate(t);
            if (originalPoint.distanceTo(reducedPoint) > tolerance) {
                throw new Error('Degree reduction exceeds tolerance; consider increasing tolerance or refining fit');
            }
        }
        return reducedCurve;
    }
}
export class NURBSSurface {
    constructor(degreeU, degreeV, controlPoints, knotsU, knotsV, weights, trimmingLoops = []) {
        this.degreeU = degreeU;
        this.degreeV = degreeV;
        this.controlPoints = controlPoints.map(row => row.map(p => p.clone()));
        this.knotsU = [...knotsU];
        this.knotsV = [...knotsV];
        this.weights = weights.map(row => [...row]);
        this.nU = controlPoints.length - 1;
        this.nV = controlPoints[0].length - 1;
        this._basisCacheU = new Map();
        this._basisCacheV = new Map();
        this._derivativeCacheU = new Map();
        this._derivativeCacheV = new Map();
        this.trimmingLoops = [];
        this._precomputedPolygons = [];
        // Array of { type, polygon: [{u, v}], bounds: {uMin, uMax, vMin, vMax} }
        this._precomputeTrimmingLoops(trimmingLoops);
    }
    evaluate(u, v) {
        const [uMin, uMax] = this.getDomainU();
        const [vMin, vMax] = this.getDomainV();
        if (u < uMin || u > uMax || v < vMin || v > vMax) {
            throw new Error(`UV parameters (${ u }, ${ v }) outside domain [${ uMin }, ${ uMax }] x [${ vMin }, ${ vMax }]`);
        }
        // Check trimming loops
        if (!this._isPointValid(u, v)) {
            return null;
        }
        // Point is outside trimmed region or inside a hole
        const spanU = this._findKnotSpan(u, this.knotsU, this.degreeU, this.nU);
        const spanV = this._findKnotSpan(v, this.knotsV, this.degreeV, this.nV);
        const basisU = this._getCachedBasis(spanU, u, this.knotsU, this.degreeU, this._basisCacheU);
        const basisV = this._getCachedBasis(spanV, v, this.knotsV, this.degreeV, this._basisCacheV);
        let x = 0, y = 0, z = 0, w = 0;
        for (let i = 0; i <= this.degreeU; i++) {
            const idxU = spanU - this.degreeU + i;
            for (let j = 0; j <= this.degreeV; j++) {
                const idxV = spanV - this.degreeV + j;
                const b = basisU[i] * basisV[j] * this.weights[idxU][idxV];
                x += this.controlPoints[idxU][idxV].x * b;
                y += this.controlPoints[idxU][idxV].y * b;
                z += this.controlPoints[idxU][idxV].z * b;
                w += b;
            }
        }
        return new Point(x / w, y / w, z / w);
    }
    derivative(u, v, orderU = 0, orderV = 0) {
        if (orderU < 0 || orderV < 0 || !Number.isInteger(orderU) || !Number.isInteger(orderV)) {
            throw new Error('Derivative orders must be non-negative integers');
        }
        if (orderU > this.degreeU || orderV > this.degreeV) {
            return new Vector(0, 0, 0);
        }
        const spanU = this._findKnotSpan(u, this.knotsU, this.degreeU, this.nU);
        const spanV = this._findKnotSpan(v, this.knotsV, this.degreeV, this.nV);
        const dersU = this._getCachedDerivatives(spanU, u, this.knotsU, this.degreeU, orderU, this._derivativeCacheU);
        const dersV = this._getCachedDerivatives(spanV, v, this.knotsV, this.degreeV, orderV, this._derivativeCacheV);
        const Aders = Array(orderU + 1).fill().map(() => Array(orderV + 1).fill().map(() => ({
            x: 0,
            y: 0,
            z: 0,
            w: 0
        })));
        for (let i = 0; i <= this.degreeU; i++) {
            const idxU = spanU - this.degreeU + i;
            for (let j = 0; j <= this.degreeV; j++) {
                const idxV = spanV - this.degreeV + j;
                const w = this.weights[idxU][idxV];
                const p = this.controlPoints[idxU][idxV];
                for (let ku = 0; ku <= orderU; ku++) {
                    for (let kv = 0; kv <= orderV; kv++) {
                        const b = dersU[ku][i] * dersV[kv][j] * w;
                        Aders[ku][kv].x += p.x * b;
                        Aders[ku][kv].y += p.y * b;
                        Aders[ku][kv].z += p.z * b;
                        Aders[ku][kv].w += b;
                    }
                }
            }
        }
        if (orderU === 0 && orderV === 0) {
            const w = Aders[0][0].w;
            return new Point(Aders[0][0].x / w, Aders[0][0].y / w, Aders[0][0].z / w);
        }
        const w = Aders[0][0].w;
        const wInv = 1 / w;
        const wInv2 = wInv * wInv;
        let result = new Vector(0, 0, 0);
        if (orderU === 1 && orderV === 0) {
            const A1 = Aders[1][0];
            const w1 = A1.w;
            result.x = (A1.x - w1 * Aders[0][0].x) * wInv2;
            result.y = (A1.y - w1 * Aders[0][0].y) * wInv2;
            result.z = (A1.z - w1 * Aders[0][0].z) * wInv2;
        } else if (orderU === 0 && orderV === 1) {
            const A1 = Aders[0][1];
            const w1 = A1.w;
            result.x = (A1.x - w1 * Aders[0][0].x) * wInv2;
            result.y = (A1.y - w1 * Aders[0][0].y) * wInv2;
            result.z = (A1.z - w1 * Aders[0][0].z) * wInv2;
        } else {
            for (let ku = 0; ku <= orderU; ku++) {
                for (let kv = 0; kv <= orderV; kv++) {
                    if (ku + kv === 0)
                        continue;
                    const coef = this._binomial(orderU, ku) * this._binomial(orderV, kv);
                    const sign = (orderU - ku + (orderV - kv)) % 2 === 0 ? 1 : -1;
                    const A = Aders[ku][kv];
                    const wTerm = Aders[orderU - ku][orderV - kv].w * coef * sign;
                    result.x += (A.x - wTerm * Aders[0][0].x) * wInv2;
                    result.y += (A.y - wTerm * Aders[0][0].y) * wInv2;
                    result.z += (A.z - wTerm * Aders[0][0].z) * wInv2;
                }
            }
        }
        return result;
    }
    _findKnotSpan(t, knots, degree, n) {
        if (t < knots[0] || t > knots[knots.length - 1]) {
            throw new Error(`Parameter t=${ t } is outside knot range [${ knots[0] }, ${ knots[knots.length - 1] }]`);
        }
        if (t === knots[knots.length - 1]) {
            return n;
        }
        let low = degree;
        let high = n + 1;
        let mid = Math.floor((low + high) / 2);
        while (t < knots[mid] || t >= knots[mid + 1]) {
            if (t < knots[mid])
                high = mid;
            else
                low = mid;
            mid = Math.floor((low + high) / 2);
        }
        return mid;
    }
    _basisFunctions(span, t, knots, degree) {
        const N = [1];
        const left = [];
        const right = [];
        for (let j = 1; j <= degree; j++) {
            left[j] = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;
            let saved = 0;
            for (let r = 0; r < j; r++) {
                const temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
        return N;
    }
    _basisDerivatives(span, t, knots, degree, n) {
        const ders = Array(n + 1).fill().map(() => Array(degree + 1).fill(0));
        const ndu = Array(degree + 1).fill().map(() => Array(degree + 1).fill(0));
        const left = Array(degree + 1).fill(0);
        const right = Array(degree + 1).fill(0);
        ndu[0][0] = 1;
        for (let j = 1; j <= degree; j++) {
            left[j] = t - knots[span + 1 - j];
            right[j] = knots[span + j] - t;
            let saved = 0;
            for (let r = 0; r < j; r++) {
                ndu[j][r] = right[r + 1] + left[j - r];
                const temp = ndu[r][j - 1] / ndu[j][r];
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }
        for (let j = 0; j <= degree; j++) {
            ders[0][j] = ndu[j][degree];
        }
        for (let r = 0; r <= degree; r++) {
            let s1 = 0, s2 = 1;
            const a = Array(degree + 1).fill().map(() => Array(degree + 1).fill(0));
            a[0][0] = 1;
            for (let k = 1; k <= n; k++) {
                let d = 0;
                const rk = r - k;
                const pk = degree - k;
                if (r >= k) {
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                    d = a[s2][0] * ndu[rk][pk];
                }
                const j1 = rk >= -1 ? 1 : -rk;
                const j2 = r - 1 <= pk ? k - 1 : degree - r;
                for (let j = j1; j <= j2; j++) {
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                    d += a[s2][j] * ndu[rk + j][pk];
                }
                if (r <= pk) {
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                    d += a[s2][k] * ndu[r][pk];
                }
                ders[k][r] = d;
                [s1, s2] = [
                    s2,
                    s1
                ];
            }
        }
        let r = degree;
        for (let k = 1; k <= n; k++) {
            for (let j = 0; j <= degree; j++) {
                ders[k][j] *= r;
            }
            r *= degree - k;
        }
        return ders;
    }
    _getCachedBasis(span, t, knots, degree, cache) {
        const key = `${ span },${ t.toFixed(6) }`;
        if (!cache.has(key)) {
            cache.set(key, this._basisFunctions(span, t, knots, degree));
        }
        return cache.get(key);
    }
    getDomainU() {
        return [
            this.knotsU[this.degreeU],
            this.knotsU[this.nU + 1]
        ];
    }
    getDomainV() {
        return [
            this.knotsV[this.degreeV],
            this.knotsV[this.nV + 1]
        ];
    }
    toString() {
        return `NURBSSurface(degreeU=${ this.degreeU }, degreeV=${ this.degreeV }, nU=${ this.nU }, nV=${ this.nV }, loops=${ this.trimmingLoops.length })`;
    }
    _getCachedDerivatives(span, t, knots, degree, order, cache) {
        const key = `${ span },${ t.toFixed(6) },${ order }`;
        if (!cache.has(key)) {
            cache.set(key, this._basisDerivatives(span, t, knots, degree, order));
        }
        return cache.get(key);
    }
    _binomial(n, k) {
        if (k < 0 || k > n)
            return 0;
        let result = 1;
        for (let i = 0; i < k; i++) {
            result *= (n - i) / (i + 1);
        }
        return result;
    }
    trim(u1, u2, v1, v2) {
        const [uMin, uMax] = this.getDomainU();
        const [vMin, vMax] = this.getDomainV();
        if (u1 < uMin || u2 > uMax || u1 >= u2 || v1 < vMin || v2 > vMax || v1 >= v2) {
            throw new Error(`Trim parameters u1=${ u1 }, u2=${ u2 }, v1=${ v1 }, v2=${ v2 } must satisfy ${ uMin } <= u1 < u2 <= ${ uMax } and ${ vMin } <= v1 < v2 <= ${ vMax }`);
        }
        let knotsU = [...this.knotsU];
        let controlPoints = this.controlPoints.map(row => row.map(p => p.clone()));
        let weights = this.weights.map(row => [...row]);
        const insertKnots = (t, knots, degree, dir) => {
            const span = this._findKnotSpan(t, knots, degree, dir === 'u' ? this.nU : this.nV);
            let multiplicity = 0;
            for (let i = span; i >= 0 && Math.abs(knots[i] - t) < 0.000001; i--) {
                multiplicity++;
            }
            return degree - multiplicity;
        };
        let insertsU1 = insertKnots(u1, knotsU, this.degreeU, 'u');
        let insertsU2 = insertKnots(u2, knotsU, this.degreeU, 'u');
        for (let i = 0; i < insertsU1; i++) {
            ({
                knots: knotsU,
                points: controlPoints,
                weights
            } = this._insertKnotU(u1, knotsU, controlPoints, weights));
        }
        for (let i = 0; i < insertsU2; i++) {
            ({
                knots: knotsU,
                points: controlPoints,
                weights
            } = this._insertKnotU(u2, knotsU, controlPoints, weights));
        }
        let knotsV = [...this.knotsV];
        let insertsV1 = insertKnots(v1, knotsV, this.degreeV, 'v');
        let insertsV2 = insertKnots(v2, knotsV, this.degreeV, 'v');
        for (let i = 0; i < insertsV1; i++) {
            ({
                knots: knotsV,
                points: controlPoints,
                weights
            } = this._insertKnotV(v1, knotsV, controlPoints, weights));
        }
        for (let i = 0; i < insertsV2; i++) {
            ({
                knots: knotsV,
                points: controlPoints,
                weights
            } = this._insertKnotV(v2, knotsV, controlPoints, weights));
        }
        let startU = 0, endU = 0, startV = 0, endV = 0;
        for (let i = 0; i < knotsU.length; i++) {
            if (knotsU[i] >= u1 - 0.000001) {
                startU = i - this.degreeU;
                break;
            }
        }
        for (let i = 0; i < knotsU.length; i++) {
            if (knotsU[i] >= u2 - 0.000001) {
                endU = i - this.degreeU;
                break;
            }
        }
        for (let i = 0; i < knotsV.length; i++) {
            if (knotsV[i] >= v1 - 0.000001) {
                startV = i - this.degreeV;
                break;
            }
        }
        for (let i = 0; i < knotsV.length; i++) {
            if (knotsV[i] >= v2 - 0.000001) {
                endV = i - this.degreeV;
                break;
            }
        }
        const trimmedKnotsU = knotsU.slice(startU, endU + this.degreeU + 1);
        const trimmedKnotsV = knotsV.slice(startV, endV + this.degreeV + 1);
        const trimmedPoints = controlPoints.slice(startU, endU + 1).map(row => row.slice(startV, endV + 1));
        const trimmedWeights = weights.slice(startU, endU + 1).map(row => row.slice(startV, endV + 1));
        return new NURBSSurface(this.degreeU, this.degreeV, trimmedPoints, trimmedKnotsU, trimmedKnotsV, trimmedWeights);
    }
    _insertKnotU(u, knotsU, controlPoints, weights) {
        const span = this._findKnotSpan(u, knotsU, this.degreeU, this.nU);
        const newKnotsU = [
            ...knotsU.slice(0, span + 1),
            u,
            ...knotsU.slice(span + 1)
        ];
        const newPoints = [];
        const newWeights = [];
        for (let i = 0; i <= span - this.degreeU; i++) {
            newPoints[i] = controlPoints[i].map(p => p.clone());
            newWeights[i] = [...weights[i]];
        }
        for (let i = span - this.degreeU + 1; i < controlPoints.length; i++) {
            newPoints[i + 1] = controlPoints[i].map(p => p.clone());
            newWeights[i + 1] = [...weights[i]];
        }
        const Q = [];
        const W = [];
        for (let i = 0; i <= this.degreeU; i++) {
            const idx = span - this.degreeU + i;
            const alpha = (u - knotsU[idx]) / (knotsU[idx + this.degreeU] - knotsU[idx]);
            Q[i] = controlPoints[idx - 1].map((p, j) => p.scale(1 - alpha).add(controlPoints[idx][j].scale(alpha)));
            W[i] = weights[idx - 1].map((w, j) => w * (1 - alpha) + weights[idx][j] * alpha);
        }
        for (let i = 0; i <= this.degreeU; i++) {
            newPoints[span - this.degreeU + i] = Q[i];
            newWeights[span - this.degreeU + i] = W[i];
        }
        return {
            knots: newKnotsU,
            points: newPoints,
            weights: newWeights
        };
    }
    _insertKnotV(v, knotsV, controlPoints, weights) {
        const span = this._findKnotSpan(v, knotsV, this.degreeV, this.nV);
        const newKnotsV = [
            ...knotsV.slice(0, span + 1),
            v,
            ...knotsV.slice(span + 1)
        ];
        const newPoints = controlPoints.map(row => row.map(p => p.clone()));
        const newWeights = weights.map(row => [...row]);
        const Q = [];
        const W = [];
        for (let j = 0; j <= this.degreeV; j++) {
            const idx = span - this.degreeV + j;
            const alpha = (v - knotsV[idx]) / (knotsV[idx + this.degreeV] - knotsV[idx]);
            Q[j] = [];
            W[j] = [];
            for (let i = 0; i < controlPoints.length; i++) {
                Q[j][i] = controlPoints[i][idx - 1].scale(1 - alpha).add(controlPoints[i][idx].scale(alpha));
                W[j][i] = weights[i][idx - 1] * (1 - alpha) + weights[i][idx] * alpha;
            }
        }
        for (let j = 0; j <= this.degreeV; j++) {
            for (let i = 0; i < controlPoints.length; i++) {
                newPoints[i][span - this.degreeV + j] = Q[j][i];
                newWeights[i][span - this.degreeV + j] = W[j][i];
            }
        }
        return {
            knots: newKnotsV,
            points: newPoints,
            weights: newWeights
        };
    }
    clearCache() {
        this._basisCacheU.clear();
        this._basisCacheV.clear();
        this._derivativeCacheU.clear();
        this._derivativeCacheV.clear();
    }
    _computeAABB(uStart, uEnd, vStart, vEnd, samples = 5) {
        const du = (uEnd - uStart) / (samples - 1);
        const dv = (vEnd - vStart) / (samples - 1);
        let minX = Infinity, minY = Infinity, minZ = Infinity;
        let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
        for (let i = 0; i < samples; i++) {
            const u = uStart + i * du;
            for (let j = 0; j < samples; j++) {
                const v = vStart + j * dv;
                const p = this.evaluate(u, v);
                minX = Math.min(minX, p.x);
                minY = Math.min(minY, p.y);
                minZ = Math.min(minZ, p.z);
                maxX = Math.max(maxX, p.x);
                maxY = Math.max(maxY, p.y);
                maxZ = Math.max(maxZ, p.z);
            }
        }
        return {
            min: new Point(minX, minY, minZ),
            max: new Point(maxX, maxY, maxZ)
        };
    }
    _aabbIntersect(aabb1, aabb2) {
        return aabb1.min.x <= aabb2.max.x && aabb1.max.x >= aabb2.min.x && aabb1.min.y <= aabb2.max.y && aabb1.max.y >= aabb2.min.y && aabb1.min.z <= aabb2.max.z && aabb1.max.z >= aabb2.min.z;
    }
    intersectSurface(otherSurface, tolerance = 0.000001, maxIterations = 100) {
        const [uMin1, uMax1] = this.getDomainU();
        const [vMin1, vMax1] = this.getDomainV();
        const [uMin2, uMax2] = otherSurface.getDomainU();
        const [vMin2, vMax2] = otherSurface.getDomainV();
        const seeds = [];
        const subdivide = (u1a, u1b, v1a, v1b, u2a, u2b, v2a, v2b, depth = 0) => {
            if (depth > 10)
                return;
            const aabb1 = this._computeAABB(u1a, u1b, v1a, v1b);
            const aabb2 = otherSurface._computeAABB(u2a, u2b, v2a, v2b);
            if (!this._aabbIntersect(aabb1, aabb2))
                return;
            const p1a = this.evaluate(u1a, v1a);
            const p1b = this.evaluate(u1b, v1b);
            const p2a = otherSurface.evaluate(u2a, v2a);
            const p2b = otherSurface.evaluate(u2b, v2b);
            const d1 = p1a.distanceTo(p1b);
            const d2 = p2a.distanceTo(p2b);
            if (d1 < tolerance && d2 < tolerance) {
                const u1 = (u1a + u1b) / 2;
                const v1 = (v1a + v1b) / 2;
                const u2 = (u2a + u2b) / 2;
                const v2 = (v2a + v2b) / 2;
                if (this.evaluate(u1, v1).distanceTo(otherSurface.evaluate(u2, v2)) < tolerance) {
                    seeds.push({
                        u1,
                        v1,
                        u2,
                        v2
                    });
                }
                return;
            }
            const u1m = (u1a + u1b) / 2;
            const v1m = (v1a + v1b) / 2;
            const u2m = (u2a + u2b) / 2;
            const v2m = (v2a + v2b) / 2;
            subdivide(u1a, u1m, v1a, v1m, u2a, u2m, v2a, v2m, depth + 1);
            subdivide(u1a, u1m, v1m, v1b, u2a, u2m, v2m, v2b, depth + 1);
            subdivide(u1m, u1b, v1a, v1m, u2m, u2b, v2a, v2m, depth + 1);
            subdivide(u1m, u1b, v1m, v1b, u2m, u2b, v2m, v2b, depth + 1);
        };
        subdivide(uMin1, uMax1, vMin1, vMax1, uMin2, uMax2, vMin2, vMax2);
        const curves = [];
        const usedSeeds = new Set();
        for (let seed of seeds) {
            const seedKey = `${ seed.u1.toFixed(6) },${ seed.v1.toFixed(6) },${ seed.u2.toFixed(6) },${ seed.v2.toFixed(6) }`;
            if (usedSeeds.has(seedKey))
                continue;
            const traced = this._traceSurfaceIntersection(seed.u1, seed.v1, seed.u2, seed.v2, otherSurface, tolerance, maxIterations);
            if (traced && traced.points.length >= 2) {
                const curve = this._fitCurve(traced.points);
                curves.push({
                    curve,
                    params1: traced.params1,
                    params2: traced.params2
                });
                // Mark seeds as used based on proximity
                for (let i = 0; i < traced.params1.length; i++) {
                    const u1 = traced.params1[i].u, v1 = traced.params1[i].v;
                    const u2 = traced.params2[i].u, v2 = traced.params2[i].v;
                    usedSeeds.add(`${ u1.toFixed(6) },${ v1.toFixed(6) },${ u2.toFixed(6) },${ v2.toFixed(6) }`);
                }
            }
        }
        return curves;
    }
    intersectCurve(curve, tolerance = 0.000001, maxIterations = 100) {
        const [tMin, tMax] = curve.getDomain();
        const [uMin, uMax] = this.getDomainU();
        const [vMin, vMax] = this.getDomainV();
        const seeds = [];
        const subdivide = (ta, tb, ua, ub, va, vb, depth = 0) => {
            if (depth > 10)
                return;
            const aabbCurve = curve._computeAABB(ta, tb);
            const aabbSurface = this._computeAABB(ua, ub, va, vb);
            if (!this._aabbIntersect(aabbCurve, aabbSurface))
                return;
            const pCa = curve.evaluate(ta);
            const pCb = curve.evaluate(tb);
            const dC = pCa.distanceTo(pCb);
            const pSa = this.evaluate(ua, va);
            const pSb = this.evaluate(ub, vb);
            const dS = pSa.distanceTo(pSb);
            if (dC < tolerance && dS < tolerance) {
                const t = (ta + tb) / 2;
                const u = (ua + ub) / 2;
                const v = (va + vb) / 2;
                if (curve.evaluate(t).distanceTo(this.evaluate(u, v)) < tolerance) {
                    seeds.push({
                        t,
                        u,
                        v
                    });
                }
                return;
            }
            const tm = (ta + tb) / 2;
            const um = (ua + ub) / 2;
            const vm = (va + vb) / 2;
            subdivide(ta, tm, ua, um, va, vm, depth + 1);
            subdivide(ta, tm, um, ub, va, vm, depth + 1);
            subdivide(ta, tm, ua, um, vm, vb, depth + 1);
            subdivide(ta, tm, um, ub, vm, vb, depth + 1);
            subdivide(tm, tb, ua, um, va, vm, depth + 1);
            subdivide(tm, tb, um, ub, va, vm, depth + 1);
            subdivide(tm, tb, ua, um, vm, vb, depth + 1);
            subdivide(tm, tb, um, ub, vm, vb, depth + 1);
        };
        subdivide(tMin, tMax, uMin, uMax, vMin, vMax);
        const curves = [];
        const usedSeeds = new Set();
        for (let seed of seeds) {
            const seedKey = `${ seed.t.toFixed(6) },${ seed.u.toFixed(6) },${ seed.v.toFixed(6) }`;
            if (usedSeeds.has(seedKey))
                continue;
            const traced = this._traceCurveIntersection(seed.t, seed.u, seed.v, curve, tolerance, maxIterations);
            if (traced && traced.points.length >= 2) {
                const curveFit = this._fitCurve(traced.points);
                curves.push({
                    curve: curveFit,
                    paramsC: traced.paramsC,
                    paramsS: traced.paramsS
                });
                for (let i = 0; i < traced.paramsC.length; i++) {
                    const t = traced.paramsC[i].t;
                    const u = traced.paramsS[i].u, v = traced.paramsS[i].v;
                    usedSeeds.add(`${ t.toFixed(6) },${ u.toFixed(6) },${ v.toFixed(6) }`);
                }
            }
        }
        return curves;
    }
    _solve3x3(matrix, b, det) {
        const [a, bRow, c] = matrix;
        const [d, e, f] = b;
        const x = (bRow[1] * (c[2] * b[0] - c[0] * b[2]) - bRow[2] * (c[1] * b[0] - c[0] * b[1]) + b[0] * (e * f - d * f)) / det;
        const y = (bRow[2] * (c[0] * b[1] - c[1] * b[0]) - b[0] * (c[2] * b[1] - c[1] * b[2]) + b[1] * (d * f - e * f)) / det;
        const z = (b[0] * (c[1] * b[2] - c[2] * b[1]) - b[1] * (c[0] * b[2] - c[2] * b[0]) + b[2] * (e * d - d * e)) / det;
        return [
            x,
            y,
            z
        ];
    }
    _solve4x4(matrix, b, det) {
        const [a, bRow, c, d] = matrix;
        const [e, f, g, h] = b;
        const x = (e * (bRow[1] * (c[2] * d[3] - c[3] * d[2]) - bRow[2] * (c[1] * d[3] - c[3] * d[1]) + bRow[3] * (c[1] * d[2] - c[2] * d[1])) - f * (bRow[0] * (c[2] * d[3] - c[3] * d[2]) - bRow[2] * (c[0] * d[3] - c[3] * d[0]) + bRow[3] * (c[0] * d[2] - c[2] * d[0])) + g * (bRow[0] * (c[1] * d[3] - c[3] * d[1]) - bRow[1] * (c[0] * d[3] - c[3] * d[0]) + bRow[3] * (c[0] * d[1] - c[1] * d[0])) - h * (bRow[0] * (c[1] * d[2] - c[2] * d[1]) - bRow[1] * (c[0] * d[2] - c[2] * d[0]) + bRow[2] * (c[0] * d[1] - c[1] * d[0]))) / det;
        const y = (f * (a[0] * (c[2] * d[3] - c[3] * d[2]) - a[2] * (c[0] * d[3] - c[3] * d[0]) + a[3] * (c[0] * d[2] - c[2] * d[0])) - e * (a[1] * (c[2] * d[3] - c[3] * d[2]) - a[2] * (c[1] * d[3] - c[3] * d[1]) + a[3] * (c[1] * d[2] - c[2] * d[1])) + g * (a[0] * (c[1] * d[3] - c[3] * d[1]) - a[1] * (c[0] * d[3] - c[3] * d[0]) + a[3] * (c[0] * d[1] - c[1] * d[0])) - h * (a[0] * (c[1] * d[2] - c[2] * d[1]) - a[1] * (c[0] * d[2] - c[2] * d[0]) + a[2] * (c[0] * d[1] - c[1] * d[0]))) / det;
        const z = (e * (a[1] * (bRow[2] * d[3] - bRow[3] * d[2]) - a[2] * (bRow[1] * d[3] - bRow[3] * d[1]) + a[3] * (bRow[1] * d[2] - bRow[2] * d[1])) - f * (a[0] * (bRow[2] * d[3] - bRow[3] * d[2]) - a[2] * (bRow[0] * d[3] - bRow[3] * d[0]) + a[3] * (bRow[0] * d[2] - bRow[2] * d[0])) + g * (a[0] * (bRow[1] * d[3] - bRow[3] * d[1]) - a[1] * (bRow[0] * d[3] - bRow[3] * d[0]) + a[3] * (bRow[0] * d[1] - bRow[1] * d[0])) - h * (a[0] * (bRow[1] * d[2] - bRow[2] * d[1]) - a[1] * (bRow[0] * d[2] - bRow[2] * d[0]) + a[2] * (bRow[0] * d[1] - bRow[1] * d[0]))) / det;
        const w = (e * (a[1] * (bRow[2] * c[3] - bRow[3] * c[2]) - a[2] * (bRow[1] * c[3] - bRow[3] * c[1]) + a[3] * (bRow[1] * c[2] - bRow[2] * c[1])) - f * (a[0] * (bRow[2] * c[3] - bRow[3] * c[2]) - a[2] * (bRow[0] * c[3] - bRow[3] * c[0]) + a[3] * (bRow[0] * c[2] - bRow[2] * c[0])) + g * (a[0] * (bRow[1] * c[3] - bRow[3] * c[1]) - a[1] * (bRow[0] * c[3] - bRow[3] * c[0]) + a[3] * (bRow[0] * c[1] - bRow[1] * c[0])) - h * (a[0] * (bRow[1] * c[2] - bRow[2] * c[1]) - a[1] * (bRow[0] * c[2] - bRow[2] * c[0]) + a[2] * (bRow[0] * c[1] - bRow[1] * c[0]))) / det;
        return [
            x,
            y,
            z,
            w
        ];
    }
    _fitCurve(points, degree = 2) {
        const n = points.length - 1;
        const knots = Array(n + degree + 2).fill(0);
        for (let i = degree + 1; i <= n; i++) {
            knots[i] = (i - degree) / (n - degree + 1);
        }
        for (let i = n + 1; i < knots.length; i++) {
            knots[i] = 1;
        }
        const weights = Array(points.length).fill(1);
        return new NURBSCurve(degree, points, knots, weights);
    }
    splitAtU(u) {
        const [uMin, uMax] = this.getDomainU();
        if (u <= uMin || u >= uMax) {
            throw new Error(`Split parameter u=${ u } must be within domain (${ uMin }, ${ uMax })`);
        }
        let knotsU = [...this.knotsU];
        let controlPoints = this.controlPoints.map(row => row.map(p => p.clone()));
        let weights = this.weights.map(row => [...row]);
        const span = this._findKnotSpan(u, knotsU, this.degreeU, this.nU);
        let multiplicity = 0;
        for (let i = span; i >= 0 && Math.abs(knotsU[i] - u) < 0.000001; i--) {
            multiplicity++;
        }
        const insertsNeeded = this.degreeU - multiplicity;
        for (let i = 0; i < insertsNeeded; i++) {
            ({
                knots: knotsU,
                points: controlPoints,
                weights
            } = this._insertKnotU(u, knotsU, controlPoints, weights));
        }
        let splitIdx = 0;
        for (let i = 0; i < knotsU.length; i++) {
            if (knotsU[i] >= u - 0.000001) {
                splitIdx = i - this.degreeU;
                break;
            }
        }
        const knotsU1 = knotsU.slice(0, splitIdx + this.degreeU + 1);
        const knotsU2 = knotsU.slice(splitIdx);
        const points1 = controlPoints.slice(0, splitIdx + 1);
        const points2 = controlPoints.slice(splitIdx);
        const weights1 = weights.slice(0, splitIdx + 1);
        const weights2 = weights.slice(splitIdx);
        return [
            new NURBSSurface(this.degreeU, this.degreeV, points1, knotsU1, [...this.knotsV], weights1),
            new NURBSSurface(this.degreeU, this.degreeV, points2, knotsU2, [...this.knotsV], weights2)
        ];
    }
    splitAtV(v) {
        const [vMin, vMax] = this.getDomainV();
        if (v <= vMin || v >= vMax) {
            throw new Error(`Split parameter v=${ v } must be within domain (${ vMin }, ${ vMax })`);
        }
        let knotsV = [...this.knotsV];
        let controlPoints = this.controlPoints.map(row => row.map(p => p.clone()));
        let weights = this.weights.map(row => [...row]);
        const span = this._findKnotSpan(v, knotsV, this.degreeV, this.nV);
        let multiplicity = 0;
        for (let i = span; i >= 0 && Math.abs(knotsV[i] - v) < 0.000001; i--) {
            multiplicity++;
        }
        const insertsNeeded = this.degreeV - multiplicity;
        for (let i = 0; i < insertsNeeded; i++) {
            ({
                knots: knotsV,
                points: controlPoints,
                weights
            } = this._insertKnotV(v, knotsV, controlPoints, weights));
        }
        let splitIdx = 0;
        for (let i = 0; i < knotsV.length; i++) {
            if (knotsV[i] >= v - 0.000001) {
                splitIdx = i - this.degreeV;
                break;
            }
        }
        const knotsV1 = knotsV.slice(0, splitIdx + this.degreeV + 1);
        const knotsV2 = knotsV.slice(splitIdx);
        const points1 = controlPoints.map(row => row.slice(0, splitIdx + 1));
        const points2 = controlPoints.map(row => row.slice(splitIdx));
        const weights1 = weights.map(row => row.slice(0, splitIdx + 1));
        const weights2 = weights.map(row => row.slice(splitIdx));
        return [
            new NURBSSurface(this.degreeU, this.degreeV, points1, [...this.knotsU], knotsV1, weights1),
            new NURBSSurface(this.degreeU, this.degreeV, points2, [...this.knotsU], knotsV2, weights2)
        ];
    }
    _refineSurfaceIntersection(u1, v1, u2, v2, otherSurface, tolerance, maxIterations) {
        const [uMin1, uMax1] = this.getDomainU();
        const [vMin1, vMax1] = this.getDomainV();
        const [uMin2, uMax2] = otherSurface.getDomainU();
        const [vMin2, vMax2] = otherSurface.getDomainV();
        let u1Refined = u1, v1Refined = v1, u2Refined = u2, v2Refined = v2;
        for (let i = 0; i < maxIterations; i++) {
            const P1 = this.evaluate(u1Refined, v1Refined);
            const P2 = otherSurface.evaluate(u2Refined, v2Refined);
            const V = P1.subtract(P2);
            const dP1u = this.derivative(u1Refined, v1Refined, 1, 0);
            const dP1v = this.derivative(u1Refined, v1Refined, 0, 1);
            const dP2u = otherSurface.derivative(u2Refined, v2Refined, 1, 0);
            const dP2v = otherSurface.derivative(u2Refined, v2Refined, 0, 1);
            const F = [
                V.x,
                V.y,
                V.z
            ];
            const J = [
                [
                    dP1u.x,
                    dP1v.x,
                    -dP2u.x,
                    -dP2v.x
                ],
                [
                    dP1u.y,
                    dP1v.y,
                    -dP2u.y,
                    -dP2v.y
                ],
                [
                    dP1u.z,
                    dP1v.z,
                    -dP2u.z,
                    -dP2v.z
                ]
            ];
            const JTJ = [
                [
                    J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0],
                    J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1],
                    J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2],
                    J[0][0] * J[0][3] + J[1][0] * J[1][3] + J[2][0] * J[2][3]
                ],
                [
                    J[0][1] * J[0][0] + J[1][1] * J[1][0] + J[2][1] * J[2][0],
                    J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1],
                    J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2],
                    J[0][1] * J[0][3] + J[1][1] * J[1][3] + J[2][1] * J[2][3]
                ],
                [
                    J[0][2] * J[0][0] + J[1][2] * J[1][0] + J[2][2] * J[2][0],
                    J[0][2] * J[0][1] + J[1][2] * J[1][1] + J[2][2] * J[2][1],
                    J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2],
                    J[0][2] * J[0][3] + J[1][2] * J[1][3] + J[2][2] * J[2][3]
                ],
                [
                    J[0][3] * J[0][0] + J[1][3] * J[1][0] + J[2][3] * J[2][0],
                    J[0][3] * J[0][1] + J[1][3] * J[1][1] + J[2][3] * J[2][1],
                    J[0][3] * J[0][2] + J[1][3] * J[1][2] + J[2][3] * J[2][2],
                    J[0][3] * J[0][3] + J[1][3] * J[1][3] + J[2][3] * J[2][3]
                ]
            ];
            const JTF = [
                J[0][0] * F[0] + J[1][0] * F[1] + J[2][0] * F[2],
                J[0][1] * F[0] + J[1][1] * F[1] + J[2][1] * F[2],
                J[0][2] * F[0] + J[1][2] * F[1] + J[2][2] * F[2],
                J[0][3] * F[0] + J[1][3] * F[1] + J[2][3] * F[2]
            ];
            const det = JTJ[0][0] * (JTJ[1][1] * JTJ[2][2] * JTJ[3][3] + JTJ[1][2] * JTJ[2][3] * JTJ[3][1] + JTJ[1][3] * JTJ[2][1] * JTJ[3][2] - JTJ[1][3] * JTJ[2][2] * JTJ[3][1] - JTJ[1][1] * JTJ[2][3] * JTJ[3][2] - JTJ[1][2] * JTJ[2][1] * JTJ[3][3]) - JTJ[0][1] * (JTJ[1][0] * JTJ[2][2] * JTJ[3][3] + JTJ[1][2] * JTJ[2][3] * JTJ[3][0] + JTJ[1][3] * JTJ[2][0] * JTJ[3][2] - JTJ[1][3] * JTJ[2][2] * JTJ[3][0] - JTJ[1][0] * JTJ[2][3] * JTJ[3][2] - JTJ[1][2] * JTJ[2][0] * JTJ[3][3]) + JTJ[0][2] * (JTJ[1][0] * JTJ[2][1] * JTJ[3][3] + JTJ[1][1] * JTJ[2][3] * JTJ[3][0] + JTJ[1][3] * JTJ[2][0] * JTJ[3][1] - JTJ[1][3] * JTJ[2][1] * JTJ[3][0] - JTJ[1][0] * JTJ[2][3] * JTJ[3][1] - JTJ[1][1] * JTJ[2][0] * JTJ[3][3]) - JTJ[0][3] * (JTJ[1][0] * JTJ[2][1] * JTJ[3][2] + JTJ[1][1] * JTJ[2][2] * JTJ[3][0] + JTJ[1][2] * JTJ[2][0] * JTJ[3][1] - JTJ[1][2] * JTJ[2][1] * JTJ[3][0] - JTJ[1][0] * JTJ[2][2] * JTJ[3][1] - JTJ[1][1] * JTJ[2][0] * JTJ[3][2]);
            if (Math.abs(det) < tolerance)
                break;
            const delta = this._solve4x4(JTJ, JTF, det);
            u1Refined -= delta[0];
            v1Refined -= delta[1];
            u2Refined -= delta[2];
            v2Refined -= delta[3];
            u1Refined = Math.max(uMin1, Math.min(uMax1, u1Refined));
            v1Refined = Math.max(vMin1, Math.min(vMax1, v1Refined));
            u2Refined = Math.max(uMin2, Math.min(uMax2, u2Refined));
            v2Refined = Math.max(vMin2, Math.min(vMax2, v2Refined));
            if (delta.every(d => Math.abs(d) < tolerance)) {
                if (this.evaluate(u1Refined, v1Refined).distanceTo(otherSurface.evaluate(u2Refined, v2Refined)) < tolerance) {
                    return {
                        u1: u1Refined,
                        v1: v1Refined,
                        u2: u2Refined,
                        v2: v2Refined
                    };
                }
                break;
            }
        }
        return null;
    }
    _refineCurveIntersection(t, u, v, curve, tolerance, maxIterations) {
        const [tMin, tMax] = curve.getDomain();
        const [uMin, uMax] = this.getDomainU();
        const [vMin, vMax] = this.getDomainV();
        let tRefined = t, uRefined = u, vRefined = v;
        for (let i = 0; i < maxIterations; i++) {
            const P1 = curve.evaluate(tRefined);
            const P2 = this.evaluate(uRefined, vRefined);
            const V = P1.subtract(P2);
            const dC = curve.derivative(tRefined, 1);
            const dSu = this.derivative(uRefined, vRefined, 1, 0);
            const dSv = this.derivative(uRefined, vRefined, 0, 1);
            const F = [
                V.x,
                V.y,
                V.z
            ];
            const J = [
                [
                    dC.x,
                    -dSu.x,
                    -dSv.x
                ],
                [
                    dC.y,
                    -dSu.y,
                    -dSv.y
                ],
                [
                    dC.z,
                    -dSu.z,
                    -dSv.z
                ]
            ];
            const JTJ = [
                [
                    J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0],
                    J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1],
                    J[0][0] * J[0][2] + J[1][0] * J[1][2] + J[2][0] * J[2][2]
                ],
                [
                    J[0][1] * J[0][0] + J[1][1] * J[1][0] + J[2][1] * J[2][0],
                    J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1],
                    J[0][1] * J[0][2] + J[1][1] * J[1][2] + J[2][1] * J[2][2]
                ],
                [
                    J[0][2] * J[0][0] + J[1][2] * J[1][0] + J[2][2] * J[2][0],
                    J[0][2] * J[0][1] + J[1][2] * J[1][1] + J[2][2] * J[2][1],
                    J[0][2] * J[0][2] + J[1][2] * J[1][2] + J[2][2] * J[2][2]
                ]
            ];
            const JTF = [
                J[0][0] * F[0] + J[1][0] * F[1] + J[2][0] * F[2],
                J[0][1] * F[0] + J[1][1] * F[1] + J[2][1] * F[2],
                J[0][2] * F[0] + J[1][2] * F[1] + J[2][2] * F[2]
            ];
            const det = JTJ[0][0] * (JTJ[1][1] * JTJ[2][2] - JTJ[1][2] * JTJ[2][1]) - JTJ[0][1] * (JTJ[1][0] * JTJ[2][2] - JTJ[1][2] * JTJ[2][0]) + JTJ[0][2] * (JTJ[1][0] * JTJ[2][1] - JTJ[1][1] * JTJ[2][0]);
            if (Math.abs(det) < tolerance)
                break;
            const delta = this._solve3x3(JTJ, JTF, det);
            tRefined -= delta[0];
            uRefined -= delta[1];
            vRefined -= delta[2];
            tRefined = Math.max(tMin, Math.min(tMax, tRefined));
            uRefined = Math.max(uMin, Math.min(uMax, uRefined));
            vRefined = Math.max(vMin, Math.min(vMax, vRefined));
            if (delta.every(d => Math.abs(d) < tolerance)) {
                if (curve.evaluate(tRefined).distanceTo(this.evaluate(uRefined, vRefined)) < tolerance) {
                    return {
                        t: tRefined,
                        u: uRefined,
                        v: vRefined
                    };
                }
                break;
            }
        }
        return null;
    }
    _traceSurfaceIntersection(seedU1, seedV1, seedU2, seedV2, otherSurface, tolerance, maxIterations, stepSize = 0.01) {
        const [uMin1, uMax1] = this.getDomainU();
        const [vMin1, vMax1] = this.getDomainV();
        const [uMin2, uMax2] = otherSurface.getDomainU();
        const [vMin2, vMax2] = otherSurface.getDomainV();
        let points = [];
        let params1 = [];
        let params2 = [];
        let u1 = seedU1, v1 = seedV1, u2 = seedU2, v2 = seedV2;
        // Refine initial seed
        let refined = this._refineSurfaceIntersection(u1, v1, u2, v2, otherSurface, tolerance, maxIterations);
        if (!refined)
            return null;
        u1 = refined.u1;
        v1 = refined.v1;
        u2 = refined.u2;
        v2 = refined.v2;
        points.push(this.evaluate(u1, v1));
        params1.push({
            u: u1,
            v: v1
        });
        params2.push({
            u: u2,
            v: v2
        });
        // Directions to trace (forward and backward)
        const directions = [
            1,
            -1
        ];
        for (let dir of directions) {
            let lastU1 = u1, lastV1 = v1, lastU2 = u2, lastV2 = v2;
            let closed = false;
            while (true) {
                // Compute tangents
                const dP1u = this.derivative(lastU1, lastV1, 1, 0);
                const dP1v = this.derivative(lastU1, lastV1, 0, 1);
                const dP2u = otherSurface.derivative(lastU2, lastV2, 1, 0);
                const dP2v = otherSurface.derivative(lastU2, lastV2, 0, 1);
                // Normal vectors (cross products)
                const n1 = dP1u.cross(dP1v);
                const n2 = dP2u.cross(dP2v);
                // Intersection direction (cross of normals)
                const dirVec = n1.cross(n2);
                if (dirVec.magnitude() < tolerance)
                    break;
                // Parallel surfaces
                const step = dirVec.normalize().scale(dir * stepSize);
                // Predict next point (simplified Euler step in 3D space, then project back)
                const nextP = points[points.length - 1].add(step);
                let nextU1 = lastU1 + stepSize * dir * (dP1u.dot(step) / dP1u.magnitude());
                let nextV1 = lastV1 + stepSize * dir * (dP1v.dot(step) / dP1v.magnitude());
                let nextU2 = lastU2 + stepSize * dir * (dP2u.dot(step) / dP2u.magnitude());
                let nextV2 = lastV2 + stepSize * dir * (dP2v.dot(step) / dP2v.magnitude());
                // Refine with Newton's method
                refined = this._refineSurfaceIntersection(nextU1, nextV1, nextU2, nextV2, otherSurface, tolerance, maxIterations);
                if (!refined)
                    break;
                nextU1 = refined.u1;
                nextV1 = refined.v1;
                nextU2 = refined.u2;
                nextV2 = refined.v2;
                const nextPoint = this.evaluate(nextU1, nextV1);
                // Check for loop closure
                if (points.length > 2 && nextPoint.distanceTo(points[0]) < tolerance * 10) {
                    closed = true;
                    break;
                }
                // Check domain boundaries
                if (nextU1 <= uMin1 || nextU1 >= uMax1 || nextV1 <= vMin1 || nextV1 >= vMax1 || nextU2 <= uMin2 || nextU2 >= uMax2 || nextV2 <= vMin2 || nextV2 >= vMax2) {
                    break;
                }
                // Check if stuck
                if (nextPoint.distanceTo(points[points.length - 1]) < tolerance)
                    break;
                points.push(nextPoint);
                params1.push({
                    u: nextU1,
                    v: nextV1
                });
                params2.push({
                    u: nextU2,
                    v: nextV2
                });
                lastU1 = nextU1;
                lastV1 = nextV1;
                lastU2 = nextU2;
                lastV2 = nextV2;
            }
            if (dir === -1 && !closed) {
                points.reverse();
                params1.reverse();
                params2.reverse();
                points = points.slice(1);
                // Remove seed point to avoid duplication
                params1 = params1.slice(1);
                params2 = params2.slice(1);
            }
        }
        if (points.length >= 2) {
            return {
                points,
                params1,
                params2
            };
        }
        return null;
    }
    _traceCurveIntersection(seedT, seedU, seedV, curve, tolerance, maxIterations, stepSize = 0.01) {
        const [tMin, tMax] = curve.getDomain();
        const [uMin, uMax] = this.getDomainU();
        const [vMin, vMax] = this.getDomainV();
        let points = [];
        let paramsC = [];
        let paramsS = [];
        let t = seedT, u = seedU, v = seedV;
        // Refine initial seed
        let refined = this._refineCurveIntersection(t, u, v, curve, tolerance, maxIterations);
        if (!refined)
            return null;
        t = refined.t;
        u = refined.u;
        v = refined.v;
        points.push(curve.evaluate(t));
        paramsC.push({ t });
        paramsS.push({
            u,
            v
        });
        const directions = [
            1,
            -1
        ];
        for (let dir of directions) {
            let lastT = t, lastU = u, lastV = v;
            let closed = false;
            while (true) {
                const dC = curve.derivative(lastT, 1);
                const dSu = this.derivative(lastU, lastV, 1, 0);
                const dSv = this.derivative(lastU, lastV, 0, 1);
                const nS = dSu.cross(dSv);
                const dirVec = dC.cross(nS);
                if (dirVec.magnitude() < tolerance)
                    break;
                const step = dirVec.normalize().scale(dir * stepSize);
                const nextP = points[points.length - 1].add(step);
                let nextT = lastT + stepSize * dir * (dC.dot(step) / dC.magnitude());
                let nextU = lastU + stepSize * dir * (dSu.dot(step) / dSu.magnitude());
                let nextV = lastV + stepSize * dir * (dSv.dot(step) / dSv.magnitude());
                refined = this._refineCurveIntersection(nextT, nextU, nextV, curve, tolerance, maxIterations);
                if (!refined)
                    break;
                nextT = refined.t;
                nextU = refined.u;
                nextV = refined.v;
                const nextPoint = curve.evaluate(nextT);
                if (points.length > 2 && nextPoint.distanceTo(points[0]) < tolerance * 10) {
                    closed = true;
                    break;
                }
                if (nextT <= tMin || nextT >= tMax || nextU <= uMin || nextU >= uMax || nextV <= vMin || nextV >= vMax) {
                    break;
                }
                if (nextPoint.distanceTo(points[points.length - 1]) < tolerance)
                    break;
                points.push(nextPoint);
                paramsC.push({ t: nextT });
                paramsS.push({
                    u: nextU,
                    v: nextV
                });
                lastT = nextT;
                lastU = nextU;
                lastV = nextV;
            }
            if (dir === -1 && !closed) {
                points.reverse();
                paramsC.reverse();
                paramsS.reverse();
                points = points.slice(1);
                paramsC = paramsC.slice(1);
                paramsS = paramsS.slice(1);
            }
        }
        if (points.length >= 2) {
            return {
                points,
                paramsC,
                paramsS
            };
        }
        return null;
    }
    trimWithLoop(uvLoopCurves, innerLoops = [], tolerance = 0.000001) {
        if (uvLoopCurves.length < 3)
            throw new Error('Outer UV loop must have at least 3 curves');
        for (let i = 0; i < uvLoopCurves.length; i++) {
            const curr = uvLoopCurves[i];
            const next = uvLoopCurves[(i + 1) % uvLoopCurves.length];
            const currEnd = curr.evaluate(curr.getDomain()[1]);
            const nextStart = next.evaluate(next.getDomain()[0]);
            if (currEnd.distanceTo(nextStart) > tolerance) {
                throw new Error('Outer UV loop curves must form a closed sequence');
            }
        }
        const trimmingLoops = [{
                type: 'outer',
                curves: uvLoopCurves
            }];
        innerLoops.forEach((innerLoop, idx) => {
            if (innerLoop.length < 3)
                throw new Error(`Inner UV loop ${ idx } must have at least 3 curves`);
            for (let i = 0; i < innerLoop.length; i++) {
                const curr = innerLoop[i];
                const next = innerLoop[(i + 1) % innerLoop.length];
                const currEnd = curr.evaluate(curr.getDomain()[1]);
                const nextStart = next.evaluate(next.getDomain()[0]);
                if (currEnd.distanceTo(nextStart) > tolerance) {
                    throw new Error(`Inner UV loop ${ idx } curves must form a closed sequence`);
                }
            }
            trimmingLoops.push({
                type: 'inner',
                curves: innerLoop
            });
        });
        const newSurface = new NURBSSurface(this.degreeU, this.degreeV, this.controlPoints.map(row => row.map(p => p.clone())), [...this.knotsU], [...this.knotsV], this.weights.map(row => [...row]), trimmingLoops);
        newSurface._precomputeTrimmingLoops(trimmingLoops, tolerance);
        return newSurface;
    }
    _isPointInLoop(u, v, loopCurves, tolerance = 0.000001) {
        // Use index-based lookup instead of reference comparison
        const loopIndex = this.trimmingLoops.findIndex(loop => loop.curves === loopCurves);
        if (loopIndex === -1 || !this._precomputedPolygons[loopIndex]) {
            throw new Error('Trimming loop not found or not precomputed');
        }
        const {polygon, bounds} = this._precomputedPolygons[loopIndex];
        if (u < bounds.uMin || u > bounds.uMax || v < bounds.vMin || v > bounds.vMax) {
            return false;
        }
        let crossings = 0;
        for (let i = 0; i < polygon.length - 1; i++) {
            const p1 = polygon[i];
            const p2 = polygon[i + 1];
            if (p1.v > v !== p2.v > v && u < p1.u + (p2.u - p1.u) * (v - p1.v) / (p2.v - p1.v)) {
                crossings++;
            }
        }
        return crossings % 2 === 1;
    }
    _isPointValid(u, v, tolerance = 0.000001) {
        if (this.trimmingLoops.length === 0)
            return true;
        let insideOuter = false;
        for (let loop of this.trimmingLoops) {
            if (loop.type === 'outer') {
                if (this._isPointInLoop(u, v, loop.curves, tolerance)) {
                    insideOuter = true;
                    break;
                }
            }
        }
        if (!insideOuter)
            return false;
        for (let loop of this.trimmingLoops) {
            if (loop.type === 'inner' && this._isPointInLoop(u, v, loop.curves, tolerance)) {
                return false;
            }
        }
        // Inside a hole
        return true;
    }
    _precomputeTrimmingLoops(trimmingLoops, tolerance = 0.000001) {
        this.trimmingLoops = trimmingLoops;
        this._precomputedPolygons = trimmingLoops.map(loop => {
            const samples = 20;
            // Fixed number of samples per curve
            const polygon = [];
            loop.curves.forEach(curve => {
                const [tMin, tMax] = curve.getDomain();
                const dt = (tMax - tMin) / (samples - 1);
                for (let i = 0; i < samples; i++) {
                    const t = tMin + i * dt;
                    const p = curve.evaluate(t);
                    polygon.push({
                        u: p.x,
                        v: p.y
                    });
                }
            });
            // Ensure closure by adjusting the last point if needed
            if (polygon.length > 0) {
                const first = polygon[0];
                const last = polygon[polygon.length - 1];
                if (Math.abs(last.u - first.u) > tolerance || Math.abs(last.v - first.v) > tolerance) {
                    polygon.push({
                        u: first.u,
                        v: first.v
                    });
                }
            }
            const uMin = Math.min(...polygon.map(p => p.u));
            const uMax = Math.max(...polygon.map(p => p.u));
            const vMin = Math.min(...polygon.map(p => p.v));
            const vMax = Math.max(...polygon.map(p => p.v));
            return {
                type: loop.type,
                polygon,
                bounds: {
                    uMin,
                    uMax,
                    vMin,
                    vMax
                }
            };
        });
    }
    elevateDegree(direction = 'u') {
        if (direction !== 'u' && direction !== 'v') {
            throw new Error('Direction must be "u" or "v"');
        }
        if (direction === 'u') {
            const newDegreeU = this.degreeU + 1;
            const newControlPoints = Array(this.nU + 2).fill().map(() => Array(this.nV + 1));
            const newWeights = Array(this.nU + 2).fill().map(() => Array(this.nV + 1));
            const newKnotsU = Array(this.nU + newDegreeU + 2);
            // Elevate each row (U-direction)
            for (let j = 0; j <= this.nV; j++) {
                const rowCurve = new NURBSCurve(this.degreeU, this.controlPoints.map(row => row[j]), this.knotsU, this.weights.map(row => row[j]));
                const elevated = rowCurve.elevateDegree();
                for (let i = 0; i <= this.nU + 1; i++) {
                    newControlPoints[i][j] = elevated.controlPoints[i];
                    newWeights[i][j] = elevated.weights[i];
                }
                if (j === 0) {
                    newKnotsU.splice(0, newKnotsU.length, ...elevated.knots);
                }
            }
            return new NURBSSurface(newDegreeU, this.degreeV, newControlPoints, newKnotsU, this.knotsV, newWeights, this.trimmingLoops);
        } else {
            const newDegreeV = this.degreeV + 1;
            const newControlPoints = Array(this.nU + 1).fill().map(() => Array(this.nV + 2));
            const newWeights = Array(this.nU + 1).fill().map(() => Array(this.nV + 2));
            const newKnotsV = Array(this.nV + newDegreeV + 2);
            // Elevate each column (V-direction)
            for (let i = 0; i <= this.nU; i++) {
                const colCurve = new NURBSCurve(this.degreeV, this.controlPoints[i], this.knotsV, this.weights[i]);
                const elevated = colCurve.elevateDegree();
                for (let j = 0; j <= this.nV + 1; j++) {
                    newControlPoints[i][j] = elevated.controlPoints[j];
                    newWeights[i][j] = elevated.weights[j];
                }
                if (i === 0) {
                    newKnotsV.splice(0, newKnotsV.length, ...elevated.knots);
                }
            }
            return new NURBSSurface(this.degreeU, newDegreeV, newControlPoints, this.knotsU, newKnotsV, newWeights, this.trimmingLoops);
        }
    }
    reduceDegree(direction = 'u', tolerance = 0.000001) {
        if (direction !== 'u' && direction !== 'v') {
            throw new Error('Direction must be "u" or "v"');
        }
        if (direction === 'u') {
            if (this.degreeU <= 1) {
                throw new Error('Cannot reduce U degree below 1');
            }
            const newDegreeU = this.degreeU - 1;
            const newControlPoints = Array(this.nU).fill().map(() => Array(this.nV + 1));
            const newWeights = Array(this.nU).fill().map(() => Array(this.nV + 1));
            const newKnotsU = Array(this.nU + newDegreeU + 2);
            // Reduce degree for each V-direction curve
            for (let j = 0; j <= this.nV; j++) {
                const rowCurve = new NURBSCurve(this.degreeU, this.controlPoints.map(row => row[j]), this.knotsU, this.weights.map(row => row[j]));
                const reduced = rowCurve.reduceDegree(tolerance);
                for (let i = 0; i < this.nU; i++) {
                    newControlPoints[i][j] = reduced.controlPoints[i];
                    newWeights[i][j] = reduced.weights[i];
                }
                if (j === 0) {
                    newKnotsU.splice(0, newKnotsU.length, ...reduced.knots);
                }
            }
            return new NURBSSurface(newDegreeU, this.degreeV, newControlPoints, newKnotsU, this.knotsV, newWeights, this.trimmingLoops);
        } else {
            if (this.degreeV <= 1) {
                throw new Error('Cannot reduce V degree below 1');
            }
            const newDegreeV = this.degreeV - 1;
            const newControlPoints = Array(this.nU + 1).fill().map(() => Array(this.nV));
            const newWeights = Array(this.nU + 1).fill().map(() => Array(this.nV));
            const newKnotsV = Array(this.nV + newDegreeV + 2);
            // Reduce degree for each U-direction curve
            for (let i = 0; i <= this.nU; i++) {
                const colCurve = new NURBSCurve(this.degreeV, this.controlPoints[i], this.knotsV, this.weights[i]);
                const reduced = colCurve.reduceDegree(tolerance);
                for (let j = 0; j < this.nV; j++) {
                    newControlPoints[i][j] = reduced.controlPoints[j];
                    newWeights[i][j] = reduced.weights[j];
                }
                if (i === 0) {
                    newKnotsV.splice(0, newKnotsV.length, ...reduced.knots);
                }
            }
            return new NURBSSurface(this.degreeU, newDegreeV, newControlPoints, this.knotsU, newKnotsV, newWeights, this.trimmingLoops);
        }
    }
    static fitPoints(pointGrid, degreeU = 3, degreeV = 3) {
        if (!Array.isArray(pointGrid) || pointGrid.length < degreeU + 1 || !pointGrid.every(row => Array.isArray(row) && row.length >= degreeV + 1)) {
            throw new Error('Point grid must be at least (degreeU + 1) x (degreeV + 1)');
        }
        const nU = pointGrid.length - 1;
        const nV = pointGrid[0].length - 1;
        const knotsU = Array(nU + degreeU + 2).fill(0);
        const knotsV = Array(nV + degreeV + 2).fill(0);
        const weights = Array(nU + 1).fill().map(() => Array(nV + 1).fill(1));
        // Generate uniform knots for U
        for (let i = degreeU + 1; i <= nU; i++) {
            knotsU[i] = (i - degreeU) / (nU - degreeU + 1);
        }
        for (let i = nU + 1; i <= nU + degreeU + 1; i++) {
            knotsU[i] = 1;
        }
        // Generate uniform knots for V
        for (let j = degreeV + 1; j <= nV; j++) {
            knotsV[j] = (j - degreeV) / (nV - degreeV + 1);
        }
        for (let j = nV + 1; j <= nV + degreeV + 1; j++) {
            knotsV[j] = 1;
        }
        // Fit curves in U direction first, then V direction (simplified global interpolation)
        const tempCurves = pointGrid.map(row => NURBSCurve.fitPoints(row, degreeV));
        const controlPoints = Array(nU + 1).fill().map(() => Array(nV + 1));
        for (let i = 0; i <= nU; i++) {
            const params = tempCurves[i].controlPoints.map((_, j) => j / nV);
            const vCurve = NURBSCurve.fitPoints(tempCurves.map(c => c.controlPoints[i]), degreeU);
            for (let j = 0; j <= nV; j++) {
                controlPoints[i][j] = vCurve.controlPoints[i] || vCurve.controlPoints[vCurve.controlPoints.length - 1];
            }
        }
        return new NURBSSurface(degreeU, degreeV, controlPoints, knotsU, knotsV, weights);
    }
}
export class Vertex {
    constructor(point) {
        this.point = point.clone();
        this.id = idCounter.getNext('vertex');
    }
    toString() {
        return `Vertex(${ this.id }, ${ this.point.toString() })`;
    }
}
export class Edge {
    constructor(curve, vertex1, vertex2, tMin = null, tMax = null) {
        this.curve = curve;
        this.vertices = [
            vertex1,
            vertex2
        ];
        this.tMin = tMin !== null ? tMin : curve.getDomain()[0];
        this.tMax = tMax !== null ? tMax : curve.getDomain()[1];
        this.id = idCounter.getNext('edge');
    }
    getStartPoint() {
        return this.curve.evaluate(this.tMin);
    }
    getEndPoint() {
        return this.curve.evaluate(this.tMax);
    }
    toString() {
        return `Edge(${ this.id }, t[${ this.tMin }, ${ this.tMax }])`;
    }
}
export class Loop {
    constructor(edges) {
        this.edges = edges;
        this.id = idCounter.getNext('loop');
        this._validate();
    }
    _validate() {
        if (this.edges.length < 1)
            throw new Error('Loop must have at least one edge');
        for (let i = 0; i < this.edges.length; i++) {
            const curr = this.edges[i];
            const next = this.edges[(i + 1) % this.edges.length];
            if (!curr.vertices[1].point.equals(next.vertices[0].point)) {
                throw new Error('Loop edges must form a closed sequence');
            }
        }
    }
    toString() {
        return `Loop(${ this.id }, ${ this.edges.length } edges)`;
    }
}
export class Face {
    constructor(surface, loops) {
        this.surface = surface;
        this.loops = loops;
        this.id = idCounter.getNext('face');
    }
    toString() {
        return `Face(${ this.id }, ${ this.loops.length } loops)`;
    }
}
export class Shell {
    constructor(faces) {
        this.faces = faces;
        this.id = idCounter.getNext('shell');
    }
    toString() {
        return `Shell(${ this.id }, ${ this.faces.length } faces)`;
    }
}
export class Solid {
    constructor(shells) {
        this.shells = shells;
        this.id = idCounter.getNext('solid');
    }
    toString() {
        return `Solid(${ this.id }, ${ this.shells.length } shells)`;
    }
}
export class SurfaceSolver {
    constructor(surface) {
        this.surface = surface;
    }
    findUV(point3D, initialU = 0.5, initialV = 0.5, tolerance = 0.000001, maxIterations = 100) {
        let u = initialU;
        let v = initialV;
        const [uMin, uMax] = this.surface.getDomainU();
        const [vMin, vMax] = this.surface.getDomainV();
        for (let i = 0; i < maxIterations; i++) {
            const P = this.surface.evaluate(u, v);
            if (!P)
                return null;
            // Outside trimmed region
            const V = P.subtract(point3D);
            const dSu = this.surface.derivative(u, v, 1, 0);
            const dSv = this.surface.derivative(u, v, 0, 1);
            const F = [
                V.x,
                V.y,
                V.z
            ];
            const J = [
                [
                    dSu.x,
                    dSv.x
                ],
                [
                    dSu.y,
                    dSv.y
                ],
                [
                    dSu.z,
                    dSv.z
                ]
            ];
            // Compute J^T J and J^T F for 2D Newton (least squares)
            const JTJ = [
                [
                    J[0][0] * J[0][0] + J[1][0] * J[1][0] + J[2][0] * J[2][0],
                    J[0][0] * J[0][1] + J[1][0] * J[1][1] + J[2][0] * J[2][1]
                ],
                [
                    J[0][1] * J[0][0] + J[1][1] * J[1][0] + J[2][1] * J[2][0],
                    J[0][1] * J[0][1] + J[1][1] * J[1][1] + J[2][1] * J[2][1]
                ]
            ];
            const JTF = [
                J[0][0] * F[0] + J[1][0] * F[1] + J[2][0] * F[2],
                J[0][1] * F[0] + J[1][1] * F[1] + J[2][1] * F[2]
            ];
            const det = JTJ[0][0] * JTJ[1][1] - JTJ[0][1] * JTJ[1][0];
            if (Math.abs(det) < tolerance)
                break;
            const deltaU = (JTJ[1][1] * JTF[0] - JTJ[0][1] * JTF[1]) / det;
            const deltaV = (-JTJ[1][0] * JTF[0] + JTJ[0][0] * JTF[1]) / det;
            u -= deltaU;
            v -= deltaV;
            u = Math.max(uMin, Math.min(uMax, u));
            v = Math.max(vMin, Math.min(vMax, v));
            if (Math.abs(deltaU) < tolerance && Math.abs(deltaV) < tolerance) {
                if (this.surface.evaluate(u, v).distanceTo(point3D) < tolerance) {
                    return new Point(u, v, 0);
                }
                break;
            }
        }
        return null;
    }
}
export class SurfaceUtils {
    edgeToUVCurve(edge, surface, tolerance = 0.000001) {
        const [tMin, tMax] = [
            edge.tMin,
            edge.tMax
        ];
        const samples = 10;
        const dt = (tMax - tMin) / (samples - 1);
        const uvPoints = [];
        const solver = new SurfaceSolver(surface);
        // Evaluate the edge and map to UV space
        let prevUV = null;
        for (let i = 0; i < samples; i++) {
            const t = tMin + i * dt;
            const point3D = edge.curve.evaluate(t);
            const initialU = prevUV ? prevUV.x : i / (samples - 1);
            const initialV = prevUV ? prevUV.y : 0;
            const uv = solver.findUV(point3D, initialU, initialV, tolerance);
            if (!uv) {
                throw new Error(`Failed to map edge point at t=${ t } to UV space`);
            }
            uvPoints.push(uv);
            prevUV = uv;
        }
        // Generate a non-decreasing knot vector
        const degree = 2;
        const n = uvPoints.length - 1;
        // Number of control points minus 1
        const knotCount = n + degree + 2;
        // Total knots required
        const knots = Array(knotCount);
        // Clamp knots at start and end
        for (let i = 0; i <= degree; i++) {
            knots[i] = 0;
            knots[knotCount - 1 - i] = 1;
        }
        // Fill interior knots with strictly increasing values
        for (let i = degree + 1; i <= n; i++) {
            knots[i] = (i - degree) / (n - degree + 1);
        }
        return new NURBSCurve(degree, uvPoints, knots, Array(uvPoints.length).fill(1));
    }
}
export class FaceSplitter {
    splitEdge(edge, tValues, tolerance = 0.000001) {
        if (!tValues.length) {
            return [edge];
        }
        tValues.sort((a, b) => a - b);
        const newEdges = [];
        let currentVertex = edge.vertices[0];
        let currentT = edge.tMin;
        for (const t of tValues) {
            if (t <= currentT + tolerance || t >= edge.tMax - tolerance) {
                continue;
            }
            const newPoint = edge.curve.evaluate(t);
            const newVertex = new Vertex(newPoint);
            newEdges.push(new Edge(edge.curve, currentVertex, newVertex, currentT, t));
            currentVertex = newVertex;
            currentT = t;
        }
        newEdges.push(new Edge(edge.curve, currentVertex, edge.vertices[1], currentT, edge.tMax));
        return newEdges;
    }
    isPointInLoop(point, loop, surface, tolerance = 0.000001) {
        let crossings = 0;
        for (let i = 0; i < loop.edges.length; i++) {
            const edge = loop.edges[i];
            const p1 = edge.getStartPoint();
            const p2 = edge.getEndPoint();
            if (p1.y > point.y !== p2.y > point.y && point.x < p1.x + (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y)) {
                crossings++;
            }
        }
        return crossings % 2 === 1;
    }
    splitFaceWithCurves(face, intersectionCurves, tolerance = 0.000001) {
        const surface = face.surface;
        const originalLoops = face.loops;
        const edgeIntersections = new Map();
        const allVertices = new Set();
        const allEdges = new Set();
        const utils = new SurfaceUtils();
        console.log('Intersection curves:', intersectionCurves);
        originalLoops.forEach(loop => loop.edges.forEach(edge => allEdges.add(edge)));
        intersectionCurves.forEach(intersection => {
            const {curve, params1} = intersection;
            let intersectionPoints = [];
            if (params1 && params1.length >= 2) {
                intersectionPoints = params1.map(param => ({
                    point: surface.evaluate(param.u, param.v),
                    tCurve: param.u
                }));
            } else {
                intersectionPoints = [
                    {
                        point: curve.evaluate(0),
                        tCurve: 0
                    },
                    {
                        point: curve.evaluate(1),
                        tCurve: 1
                    }
                ];
            }
            console.log('Intersection points:', intersectionPoints);
            const uvCurve = utils.edgeToUVCurve(new Edge(curve, new Vertex(intersectionPoints[0].point), new Vertex(intersectionPoints[intersectionPoints.length - 1].point), 0, 1), surface, tolerance);
            let foundIntersections = false;
            allEdges.forEach(edge => {
                const edgeUVCurve = utils.edgeToUVCurve(edge, surface, tolerance);
                const uvIntersections = edgeUVCurve.intersect(uvCurve, tolerance);
                const tValues = uvIntersections.filter(i => i.t1 >= edge.tMin && i.t1 <= edge.tMax).map(i => i.t1);
                console.log(`UV intersections with edge ${ edge.tMin }-${ edge.tMax }:`, tValues);
                if (tValues.length > 0) {
                    edgeIntersections.set(edge, tValues.map(t => ({ tEdge: t })));
                    foundIntersections = true;
                }
            });
            let startVertex, endVertex;
            if (!foundIntersections) {
                console.log('No UV intersections, forcing split with endpoints');
                let intersectionsFound = 0;
                intersectionPoints.forEach(({point}) => {
                    allEdges.forEach(edge => {
                        const edgeCurve = edge.curve;
                        const edgeStart = edgeCurve.evaluate(edge.tMin);
                        const edgeEnd = edgeCurve.evaluate(edge.tMax);
                        if (!edgeStart || !edgeEnd)
                            return;
                        const edgeVec = edgeEnd.subtract(edgeStart);
                        if (!edgeVec || typeof edgeVec.magnitude !== 'function')
                            return;
                        const edgeLength = edgeVec.magnitude();
                        if (edgeLength < tolerance)
                            return;
                        const vec = point.subtract(edgeStart);
                        const t = vec.dot(edgeVec) / (edgeLength * edgeLength);
                        const tEdge = t * (edge.tMax - edge.tMin) + edge.tMin;
                        if (tEdge >= edge.tMin + tolerance && tEdge <= edge.tMax - tolerance) {
                            const projectedPoint = edgeCurve.evaluate(tEdge);
                            if (projectedPoint.distanceTo(point) < tolerance * 10) {
                                edgeIntersections.set(edge, (edgeIntersections.get(edge) || []).concat([{ tEdge }]));
                                intersectionsFound++;
                            }
                        }
                    });
                });
                if (intersectionsFound < 2 && intersectionPoints.length >= 2) {
                    console.log('Forcing two intersections on distinct edges');
                    const edgesArray = [...allEdges];
                    if (edgesArray.length >= 2) {
                        const edge1 = edgesArray[0];
                        const edge2 = edgesArray[2] || edgesArray[1];
                        // Prefer opposite edge (e.g., top/bottom)
                        const t1 = (edge1.tMin + edge1.tMax) / 2;
                        const t2 = (edge2.tMin + edge2.tMax) / 2;
                        edgeIntersections.set(edge1, [{ tEdge: t1 }]);
                        edgeIntersections.set(edge2, [{ tEdge: t2 }]);
                        startVertex = new Vertex(edge1.curve.evaluate(t1));
                        endVertex = new Vertex(edge2.curve.evaluate(t2));
                        allVertices.add(startVertex);
                        allVertices.add(endVertex);
                    }
                }
            }
            if (!startVertex || !endVertex) {
                startVertex = [...allVertices].find(v => v.point.distanceTo(intersectionPoints[0].point) < tolerance) || new Vertex(intersectionPoints[0].point);
                endVertex = [...allVertices].find(v => v.point.distanceTo(intersectionPoints[intersectionPoints.length - 1].point) < tolerance) || new Vertex(intersectionPoints[intersectionPoints.length - 1].point);
                allVertices.add(startVertex);
                allVertices.add(endVertex);
            }
            const newEdge = new Edge(curve, startVertex, endVertex, 0, 1);
            allEdges.add(newEdge);
        });
        edgeIntersections.forEach((intersections, edge) => {
            const tValues = intersections.map(i => i.tEdge).filter(t => t > edge.tMin + tolerance && t < edge.tMax - tolerance);
            console.log(`Splitting edge ${ edge.tMin }-${ edge.tMax } at:`, tValues);
            const newEdges = this.splitEdge(edge, tValues, tolerance);
            newEdges.forEach(e => {
                allVertices.add(e.vertices[0]);
                allVertices.add(e.vertices[1]);
                allEdges.add(e);
            });
            allEdges.delete(edge);
        });
        const loops = [];
        const visitedEdges = new Set();
        while (allEdges.size > 0) {
            const startEdge = [...allEdges][0];
            if (visitedEdges.has(startEdge)) {
                allEdges.delete(startEdge);
                continue;
            }
            const loopEdges = [];
            let currentEdge = startEdge;
            let currentVertex = currentEdge.vertices[0];
            let closed = false;
            do {
                loopEdges.push(currentEdge);
                visitedEdges.add(currentEdge);
                allEdges.delete(currentEdge);
                const nextVertex = currentEdge.vertices[0].point.equals(currentVertex.point, tolerance) ? currentEdge.vertices[1] : currentEdge.vertices[0];
                let nextEdge = null;
                for (const edge of allEdges) {
                    if (!visitedEdges.has(edge) && (edge.vertices[0].point.equals(nextVertex.point, tolerance) || edge.vertices[1].point.equals(nextVertex.point, tolerance))) {
                        nextEdge = edge;
                        break;
                    }
                }
                if (!nextEdge) {
                    if (nextVertex.point.equals(loopEdges[0].vertices[0].point, tolerance)) {
                        closed = true;
                    }
                    break;
                }
                currentEdge = nextEdge;
                currentVertex = nextVertex;
            } while (!closed && loopEdges.length < 1000);
            if (closed && loopEdges.length >= 3) {
                loops.push(new Loop(loopEdges));
            }
        }
        console.log('Loops formed:', loops.length);
        const faces = [];
        const outerLoops = [];
        const innerLoopsMap = new Map();
        loops.forEach(loop => {
            const centroid = loop.edges.reduce((sum, edge) => sum.add(edge.getStartPoint()), new Point(0, 0, 0)).scale(1 / loop.edges.length);
            let isOuter = true;
            for (const otherLoop of loops) {
                if (loop !== otherLoop && this.isPointInLoop(centroid, otherLoop, surface, tolerance)) {
                    isOuter = false;
                    break;
                }
            }
            if (isOuter) {
                outerLoops.push(loop);
                innerLoopsMap.set(loop, []);
            } else {
                for (const outer of outerLoops) {
                    if (this.isPointInLoop(centroid, outer, surface, tolerance)) {
                        innerLoopsMap.get(outer).push(loop);
                        break;
                    }
                }
            }
        });
        console.log('Outer loops:', outerLoops.length);
        outerLoops.forEach(outerLoop => {
            const innerLoops = innerLoopsMap.get(outerLoop);
            const outerUVCurves = outerLoop.edges.map(edge => utils.edgeToUVCurve(edge, surface, tolerance));
            const innerUVCurves = innerLoops.length > 0 ? innerLoops.map(loop => loop.edges.map(edge => utils.edgeToUVCurve(edge, surface, tolerance))) : [];
            const trimmedSurface = surface.trimWithLoop(outerUVCurves, innerUVCurves, tolerance);
            const faceLoops = [
                new Loop(outerLoop.edges),
                ...innerLoops.map(loop => new Loop(loop.edges))
            ];
            faces.push(new Face(trimmedSurface, faceLoops));
        });
        if (faces.length === 0) {
            throw new Error('No valid faces produced from splitting');
        }
        console.log('Faces produced:', faces.length);
        return faces;
    }
}
export class Tester {
    constructor(testsInstance) {
        this.tests = testsInstance;
    }
    runTests() {
        const results = [];
        const testMethods = Object.getOwnPropertyNames(Object.getPrototypeOf(this.tests)).filter(name => name !== 'constructor' && typeof this.tests[name] === 'function');
        for (const methodName of testMethods) {
            try {
                this.tests[methodName]();
                results.push({
                    test: methodName,
                    passed: true,
                    message: 'Test passed'
                });
            } catch (error) {
                results.push({
                    test: methodName,
                    passed: false,
                    message: `Test failed: ${ error.message }`,
                    error: error
                });
            }
        }
        this._reportResults(results);
        return results;
    }
    _reportResults(results) {
        console.log('Test Results:');
        results.forEach(result => {
            const status = result.passed ? 'PASSED' : 'FAILED';
            console.log(`${ status } - ${ result.test }: ${ result.message }`);
            if (!result.passed) {
                console.log(result.error);
            }
        });
        const total = results.length;
        const passed = results.filter(r => r.passed).length;
        console.log(`Summary: ${ passed }/${ total } tests passed`);
    }
}
export class Tests {
    testNURBSCurveEvaluation() {
        const curve = new NURBSCurve(1, [
            new Point(0, 0, 0),
            new Point(1, 1, 1)
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]);
        const point = curve.evaluate(0.5);
        if (!point || point.distanceTo(new Point(0.5, 0.5, 0.5)) > 0.000001) {
            throw new Error('NURBSCurve evaluation at t=0.5 should be (0.5, 0.5, 0.5)');
        }
    }
    testVertexCreation() {
        const vertex = new Vertex(new Point(1, 2, 3));
        if (!vertex.point.equals(new Point(1, 2, 3))) {
            throw new Error('Vertex should store the correct point');
        }
    }
    testEdgeEndPoints() {
        const v1 = new Vertex(new Point(0, 0, 0));
        const v2 = new Vertex(new Point(1, 0, 0));
        const curve = new NURBSCurve(1, [
            v1.point,
            v2.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]);
        const edge = new Edge(curve, v1, v2);
        const start = edge.getStartPoint();
        const end = edge.getEndPoint();
        if (!start.equals(v1.point) || !end.equals(v2.point)) {
            throw new Error('Edge endpoints should match vertex points');
        }
    }
    testLoopValidation() {
        const v1 = new Vertex(new Point(0, 0, 0));
        const v2 = new Vertex(new Point(1, 0, 0));
        const v3 = new Vertex(new Point(1, 1, 0));
        const e1 = new Edge(new NURBSCurve(1, [
            v1.point,
            v2.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v1, v2);
        const e2 = new Edge(new NURBSCurve(1, [
            v2.point,
            v3.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v2, v3);
        try {
            new Loop([
                e1,
                e2
            ]);
            // Not closed
            throw new Error('Loop should throw error for non-closed sequence');
        } catch (error) {
            if (error.message !== 'Loop edges must form a closed sequence') {
                throw new Error('Loop validation should detect non-closed sequence');
            }
        }
    }
    testShellConstruction() {
        const surface = new NURBSSurface(1, 1, [
            [
                new Point(0, 0, 0),
                new Point(0, 1, 0)
            ],
            [
                new Point(1, 0, 0),
                new Point(1, 1, 0)
            ]
        ], [
            0,
            0,
            1,
            1
        ], [
            0,
            0,
            1,
            1
        ], [
            [
                1,
                1
            ],
            [
                1,
                1
            ]
        ]);
        const v1 = new Vertex(new Point(0, 0, 0));
        const v2 = new Vertex(new Point(1, 0, 0));
        const v3 = new Vertex(new Point(1, 1, 0));
        const v4 = new Vertex(new Point(0, 1, 0));
        const e1 = new Edge(new NURBSCurve(1, [
            v1.point,
            v2.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v1, v2);
        const e2 = new Edge(new NURBSCurve(1, [
            v2.point,
            v3.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v2, v3);
        const e3 = new Edge(new NURBSCurve(1, [
            v3.point,
            v4.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v3, v4);
        const e4 = new Edge(new NURBSCurve(1, [
            v4.point,
            v1.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v4, v1);
        const loop = new Loop([
            e1,
            e2,
            e3,
            e4
        ]);
        const face = new Face(surface, [loop]);
        const shell = new Shell([face]);
        if (shell.faces.length !== 1 || shell.faces[0] !== face) {
            throw new Error('Shell should contain the correct face');
        }
    }
    testSolidConstruction() {
        const surface = new NURBSSurface(1, 1, [
            [
                new Point(0, 0, 0),
                new Point(0, 1, 0)
            ],
            [
                new Point(1, 0, 0),
                new Point(1, 1, 0)
            ]
        ], [
            0,
            0,
            1,
            1
        ], [
            0,
            0,
            1,
            1
        ], [
            [
                1,
                1
            ],
            [
                1,
                1
            ]
        ]);
        const v1 = new Vertex(new Point(0, 0, 0));
        const v2 = new Vertex(new Point(1, 0, 0));
        const v3 = new Vertex(new Point(1, 1, 0));
        const v4 = new Vertex(new Point(0, 1, 0));
        const e1 = new Edge(new NURBSCurve(1, [
            v1.point,
            v2.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v1, v2);
        const e2 = new Edge(new NURBSCurve(1, [
            v2.point,
            v3.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v2, v3);
        const e3 = new Edge(new NURBSCurve(1, [
            v3.point,
            v4.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v3, v4);
        const e4 = new Edge(new NURBSCurve(1, [
            v4.point,
            v1.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v4, v1);
        const loop = new Loop([
            e1,
            e2,
            e3,
            e4
        ]);
        const face = new Face(surface, [loop]);
        const shell = new Shell([face]);
        const solid = new Solid([shell]);
        if (solid.shells.length !== 1 || solid.shells[0] !== shell) {
            throw new Error('Solid should contain the correct shell');
        }
    }
    testSurfaceSolverAccuracy() {
        const surface = new NURBSSurface(1, 1, [
            [
                new Point(0, 0, 0),
                new Point(0, 1, 0)
            ],
            [
                new Point(1, 0, 0),
                new Point(1, 1, 0)
            ]
        ], [
            0,
            0,
            1,
            1
        ], [
            0,
            0,
            1,
            1
        ], [
            [
                1,
                1
            ],
            [
                1,
                1
            ]
        ]);
        const solver = new SurfaceSolver(surface);
        const target = new Point(0.5, 0.5, 0);
        const uv = solver.findUV(target);
        if (!uv || Math.abs(uv.x - 0.5) > 0.000001 || Math.abs(uv.y - 0.5) > 0.000001) {
            throw new Error('SurfaceSolver should find accurate UV coordinates');
        }
    }
    testSurfaceUtilsUVMapping() {
        const surface = new NURBSSurface(1, 1, [
            [
                new Point(0, 0, 0),
                new Point(0, 1, 0)
            ],
            [
                new Point(1, 0, 0),
                new Point(1, 1, 0)
            ]
        ], [
            0,
            0,
            1,
            1
        ], [
            0,
            0,
            1,
            1
        ], [
            [
                1,
                1
            ],
            [
                1,
                1
            ]
        ]);
        const v1 = new Vertex(new Point(0, 0, 0));
        const v2 = new Vertex(new Point(1, 0, 0));
        const edge = new Edge(new NURBSCurve(1, [
            v1.point,
            v2.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v1, v2);
        const utils = new SurfaceUtils();
        const uvCurve = utils.edgeToUVCurve(edge, surface);
        const midPoint = uvCurve.evaluate(0.5);
        if (!midPoint || Math.abs(midPoint.x - 0.5) > 0.000001 || Math.abs(midPoint.y - 0) > 0.000001) {
            throw new Error('SurfaceUtils should map edge to correct UV curve');
        }
    }
    testSurfaceEvaluationInsideTrim() {
        const surface = new NURBSSurface(1, 1, [
            [
                new Point(0, 0, 0),
                new Point(0, 1, 0)
            ],
            [
                new Point(1, 0, 0),
                new Point(1, 1, 0)
            ]
        ], [
            0,
            0,
            1,
            1
        ], [
            0,
            0,
            1,
            1
        ], [
            [
                1,
                1
            ],
            [
                1,
                1
            ]
        ]);
        const uvLoop = [
            new NURBSCurve(1, [
                new Point(0.25, 0.25, 0),
                new Point(0.75, 0.25, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ]),
            new NURBSCurve(1, [
                new Point(0.75, 0.25, 0),
                new Point(0.75, 0.75, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ]),
            new NURBSCurve(1, [
                new Point(0.75, 0.75, 0),
                new Point(0.25, 0.75, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ]),
            new NURBSCurve(1, [
                new Point(0.25, 0.75, 0),
                new Point(0.25, 0.25, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ])
        ];
        const trimmedSurface = surface.trimWithLoop(uvLoop);
        const point = trimmedSurface.evaluate(0.5, 0.5);
        if (!point || point.distanceTo(new Point(0.5, 0.5, 0)) > 0.000001) {
            throw new Error('Point inside trim should evaluate correctly');
        }
    }
    testSurfaceEvaluationOutsideTrim() {
        const surface = new NURBSSurface(1, 1, [
            [
                new Point(0, 0, 0),
                new Point(0, 1, 0)
            ],
            [
                new Point(1, 0, 0),
                new Point(1, 1, 0)
            ]
        ], [
            0,
            0,
            1,
            1
        ], [
            0,
            0,
            1,
            1
        ], [
            [
                1,
                1
            ],
            [
                1,
                1
            ]
        ]);
        const uvLoop = [
            new NURBSCurve(1, [
                new Point(0.25, 0.25, 0),
                new Point(0.75, 0.25, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ]),
            new NURBSCurve(1, [
                new Point(0.75, 0.25, 0),
                new Point(0.75, 0.75, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ]),
            new NURBSCurve(1, [
                new Point(0.75, 0.75, 0),
                new Point(0.25, 0.75, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ]),
            new NURBSCurve(1, [
                new Point(0.25, 0.75, 0),
                new Point(0.25, 0.25, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ])
        ];
        const trimmedSurface = surface.trimWithLoop(uvLoop);
        const point = trimmedSurface.evaluate(0.1, 0.1);
        if (point !== null) {
            throw new Error('Point outside trim should return null');
        }
    }
    testFaceSplitting() {
        const surface = new NURBSSurface(1, 1, [
            [
                new Point(0, 0, 0),
                new Point(0, 1, 0)
            ],
            [
                new Point(1, 0, 0),
                new Point(1, 1, 0)
            ]
        ], [
            0,
            0,
            1,
            1
        ], [
            0,
            0,
            1,
            1
        ], [
            [
                1,
                1
            ],
            [
                1,
                1
            ]
        ]);
        const v1 = new Vertex(new Point(0, 0, 0));
        const v2 = new Vertex(new Point(1, 0, 0));
        const v3 = new Vertex(new Point(1, 1, 0));
        const v4 = new Vertex(new Point(0, 1, 0));
        const e1 = new Edge(new NURBSCurve(1, [
            v1.point,
            v2.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v1, v2);
        const e2 = new Edge(new NURBSCurve(1, [
            v2.point,
            v3.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v2, v3);
        const e3 = new Edge(new NURBSCurve(1, [
            v3.point,
            v4.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v3, v4);
        const e4 = new Edge(new NURBSCurve(1, [
            v4.point,
            v1.point
        ], [
            0,
            0,
            1,
            1
        ], [
            1,
            1
        ]), v4, v1);
        const loop = new Loop([
            e1,
            e2,
            e3,
            e4
        ]);
        const face = new Face(surface, [loop]);
        // Adjusted intersection curve to ensure it intersects e1 and e3
        const intersectionCurve = {
            curve: new NURBSCurve(1, [
                new Point(0.5, 0, 0),
                new Point(0.5, 1, 0)
            ], [
                0,
                0,
                1,
                1
            ], [
                1,
                1
            ]),
            params1: [
                {
                    u: 0.5,
                    v: 0
                },
                {
                    u: 0.5,
                    v: 1
                }
            ]
        };
        const splitter = new FaceSplitter();
        const newFaces = splitter.splitFaceWithCurves(face, [intersectionCurve]);
        // Debug: Log intersections
        const intersections = e1.curve.intersect(intersectionCurve.curve);
        console.log('Intersections with e1:', intersections);
        if (newFaces.length !== 2) {
            throw new Error(`Face splitting should produce exactly 2 faces, got ${ newFaces.length }`);
        }
    }
}
export class FaceSplitterHelper {
    findEdgeCurveIntersections(edge, intersectionCurve, tolerance = 0.000001) {
        const curve = intersectionCurve.curve;
        const intersections = curve.intersect(edge.curve, tolerance) || [];
        // Ensure array
        const results = [];
        for (let intersection of intersections) {
            const tEdge = intersection.t2;
            const tCurve = intersection.t1;
            if (tEdge >= edge.tMin - tolerance && tEdge <= edge.tMax + tolerance && tCurve >= curve.tMin - tolerance && tCurve <= curve.tMax + tolerance) {
                const point = edge.curve.evaluate(tEdge);
                results.push({
                    point,
                    tEdge,
                    tCurve
                });
            }
        }
        return results;
    }
}
const idCounter = {
    vertex: 0,
    edge: 0,
    loop: 0,
    face: 0,
    shell: 0,
    solid: 0,
    getNext(type) {
        return this[type]++;
    }
};
// Function to split a face along an intersection curve
export function splitFaceWithCurve(face, intersectionCurve, tolerance = 0.000001) {
    const {curve, params1} = intersectionCurve;
    // From intersectSurface or intersectCurve
    const surface = face.surface;
    // Step 1: Create vertices along the intersection curve
    const vertices = params1.map(param => new Vertex(surface.evaluate(param.u, param.v)));
    if (vertices.length < 2)
        throw new Error('Intersection curve must have at least two points');
    // Step 2: Create edges along the intersection curve
    const edges = [];
    for (let i = 0; i < vertices.length - 1; i++) {
        const tMin = i / (vertices.length - 1);
        const tMax = (i + 1) / (vertices.length - 1);
        edges.push(new Edge(curve, vertices[i], vertices[i + 1], tMin, tMax));
    }
    // Step 3: Split the original loop (assuming a simple rectangular face for now)
    const originalLoop = face.loops[0];
    const originalEdges = originalLoop.edges;
    // Find where the intersection curve intersects the boundary (simplified assumption: intersects two edges)
    const intersectionPoints = [
        vertices[0].point,
        // Start of intersection curve
        vertices[vertices.length - 1].point
    ];
    // End of intersection curve
    let splitEdges1 = [], splitEdges2 = [];
    let addedStart = false, addedEnd = false;
    for (let edge of originalEdges) {
        const start = edge.getStartPoint();
        const end = edge.getEndPoint();
        if (start.distanceTo(intersectionPoints[0]) < tolerance && !addedStart) {
            splitEdges1.push(edge);
            splitEdges2.push(new Edge(edge.curve, edge.vertices[0], vertices[0], edge.tMin, edge.tMin));
            splitEdges1.push(edges[0]);
            // First intersection edge
            addedStart = true;
        } else if (end.distanceTo(intersectionPoints[1]) < tolerance && !addedEnd) {
            splitEdges1.push(new Edge(edge.curve, vertices[vertices.length - 1], edge.vertices[1], edge.tMax, edge.tMax));
            splitEdges2.push(edge);
            splitEdges2 = splitEdges2.concat(edges.slice().reverse());
            // Reverse intersection edges for second loop
            addedEnd = true;
        } else {
            splitEdges1.push(edge);
            splitEdges2.push(edge);
        }
    }
    // Step 4: Create new loops
    const loop1 = new Loop(splitEdges1);
    const loop2 = new Loop(splitEdges2);
    // Step 5: Trim the surface and create new faces
    const [uMin, uMax] = surface.getDomainU();
    const [vMin, vMax] = surface.getDomainV();
    const midU = (params1[0].u + params1[params1.length - 1].u) / 2;
    // Simplified split point
    const surface1 = surface.trim(uMin, midU, vMin, vMax);
    const surface2 = surface.trim(midU, uMax, vMin, vMax);
    const face1 = new Face(surface1, [loop1]);
    const face2 = new Face(surface2, [loop2]);
    return [
        face1,
        face2
    ];
}
// Helper: Find intersection points between an edge and an intersection curve
export function findEdgeCurveIntersections(edge, intersectionCurve, tolerance = 0.000001) {
    const curve = intersectionCurve.curve;
    const intersections = curve.intersectCurve(edge.curve, tolerance);
    const results = [];
    for (let intersection of intersections) {
        const tEdge = intersection.t2;
        // Parameter on edge's curve
        const tCurve = intersection.t1;
        // Parameter on intersection curve
        if (tEdge >= edge.tMin - tolerance && tEdge <= edge.tMax + tolerance && tCurve >= curve.tMin - tolerance && tCurve <= curve.tMax + tolerance) {
            const point = edge.curve.evaluate(tEdge);
            results.push({
                point,
                tEdge,
                tCurve
            });
        }
    }
    return results;
}
// Helper: Split an edge at given parameter values
export function splitEdge(edge, tValues, tolerance = 0.000001) {
    if (tValues.length === 0)
        return [edge];
    tValues.sort((a, b) => a - b);
    const newEdges = [];
    let currentVertex = edge.vertices[0];
    let currentT = edge.tMin;
    for (let t of tValues) {
        if (t <= currentT + tolerance || t >= edge.tMax - tolerance)
            continue;
        const newPoint = edge.curve.evaluate(t);
        const newVertex = new Vertex(newPoint);
        newEdges.push(new Edge(edge.curve, currentVertex, newVertex, currentT, t));
        currentVertex = newVertex;
        currentT = t;
    }
    newEdges.push(new Edge(edge.curve, currentVertex, edge.vertices[1], currentT, edge.tMax));
    return newEdges;
}
// Helper: Check if a point is inside a loop (2D projection for simplicity)
export function isPointInLoop(point, loop, surface, tolerance = 0.000001) {
    // Project to UV space for simplicity (assumes surface is roughly planar in this example)
    let crossings = 0;
    for (let i = 0; i < loop.edges.length; i++) {
        const edge = loop.edges[i];
        const p1 = edge.getStartPoint();
        const p2 = edge.getEndPoint();
        if (p1.y > point.y !== p2.y > point.y && point.x < p1.x + (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y)) {
            crossings++;
        }
    }
    return crossings % 2 === 1;
}
// Helper: Convert a 3D edge to a UV NURBS curve on the surface
export function edgeToUVCurve(edge, surface, tolerance = 0.000001) {
    const [tMin, tMax] = [
        edge.tMin,
        edge.tMax
    ];
    const samples = 10;
    const dt = (tMax - tMin) / (samples - 1);
    const uvPoints = [];
    for (let i = 0; i < samples; i++) {
        const t = tMin + i * dt;
        const point3D = edge.curve.evaluate(t);
        let minDist = Infinity;
        let bestUV = null;
        const uRange = surface.getDomainU();
        const vRange = surface.getDomainV();
        const uSteps = 20, vSteps = 20;
        for (let u = uRange[0]; u <= uRange[1]; u += (uRange[1] - uRange[0]) / uSteps) {
            for (let v = vRange[0]; v <= vRange[1]; v += (vRange[1] - vRange[0]) / vSteps) {
                const p = surface.evaluate(u, v);
                if (p && p.distanceTo(point3D) < minDist) {
                    minDist = p.distanceTo(point3D);
                    bestUV = new Point(u, v, 0);
                }
            }
        }
        if (minDist > tolerance)
            throw new Error('Could not map edge point to surface UV space accurately');
        uvPoints.push(bestUV);
    }
    return new NURBSCurve(2, uvPoints, [
        0,
        0,
        0,
        ...Array(uvPoints.length - 2).fill().map((_, i) => (i + 1) / (uvPoints.length - 1)),
        1,
        1,
        1
    ], Array(uvPoints.length).fill(1));
}
// Updated splitFaceWithCurves
export function splitFaceWithCurves(face, intersectionCurves, tolerance = 0.000001) {
    const surface = face.surface;
    const originalLoops = face.loops;
    const edgeIntersections = new Map();
    const allVertices = new Set();
    const allEdges = new Set();
    originalLoops.forEach(loop => loop.edges.forEach(edge => allEdges.add(edge)));
    intersectionCurves.forEach(intersection => {
        const curve = intersection.curve;
        allEdges.forEach(edge => {
            const intersections = findEdgeCurveIntersections(edge, intersection, tolerance);
            if (intersections.length > 0) {
                edgeIntersections.set(edge, (edgeIntersections.get(edge) || []).concat(intersections));
            }
        });
    });
    const newEdgesMap = new Map();
    edgeIntersections.forEach((intersections, edge) => {
        const tValues = intersections.map(i => i.tEdge);
        const newEdges = splitEdge(edge, tValues, tolerance);
        newEdgesMap.set(edge, newEdges);
        newEdges.forEach(e => {
            allVertices.add(e.vertices[0]);
            allVertices.add(e.vertices[1]);
            allEdges.add(e);
        });
        allEdges.delete(edge);
    });
    allEdges.forEach(edge => {
        if (!newEdgesMap.has(edge)) {
            allVertices.add(edge.vertices[0]);
            allVertices.add(edge.vertices[1]);
        }
    });
    const curveEdges = [];
    intersectionCurves.forEach(intersection => {
        const {curve, params1} = intersection;
        const vertices = params1.map(param => new Vertex(surface.evaluate(param.u, param.v)));
        for (let i = 0; i < vertices.length - 1; i++) {
            const tMin = i / (vertices.length - 1);
            const tMax = (i + 1) / (vertices.length - 1);
            const newEdge = new Edge(curve, vertices[i], vertices[i + 1], tMin, tMax);
            curveEdges.push(newEdge);
            allEdges.add(newEdge);
            allVertices.add(vertices[i]);
            allVertices.add(vertices[i + 1]);
        }
    });
    const loops = [];
    const visitedEdges = new Set();
    while (allEdges.size > 0) {
        const startEdge = [...allEdges][0];
        if (visitedEdges.has(startEdge)) {
            allEdges.delete(startEdge);
            continue;
        }
        const loopEdges = [];
        let currentEdge = startEdge;
        let currentVertex = currentEdge.vertices[0];
        let closed = false;
        do {
            loopEdges.push(currentEdge);
            visitedEdges.add(currentEdge);
            allEdges.delete(currentEdge);
            const nextVertex = currentEdge.vertices[0].point.equals(currentVertex.point) ? currentEdge.vertices[1] : currentEdge.vertices[0];
            let nextEdge = null;
            for (let edge of allEdges) {
                if (!visitedEdges.has(edge) && (edge.vertices[0].point.equals(nextVertex.point) || edge.vertices[1].point.equals(nextVertex.point))) {
                    nextEdge = edge;
                    break;
                }
            }
            if (!nextEdge) {
                if (nextVertex.point.equals(loopEdges[0].vertices[0].point)) {
                    closed = true;
                    break;
                }
                break;
            }
            currentEdge = nextEdge;
            currentVertex = nextVertex;
        } while (!closed && loopEdges.length < 1000);
        if (closed && loopEdges.length >= 3) {
            loops.push(new Loop(loopEdges));
        }
    }
    const faces = [];
    const outerLoops = [];
    const innerLoopsMap = new Map();
    loops.forEach(loop => {
        const centroid = loop.edges.reduce((sum, edge) => sum.add(edge.getStartPoint()), new Point(0, 0, 0)).scale(1 / loop.edges.length);
        let isOuter = true;
        for (let otherLoop of loops) {
            if (loop !== otherLoop && isPointInLoop(centroid, otherLoop, surface, tolerance)) {
                isOuter = false;
                break;
            }
        }
        if (isOuter) {
            outerLoops.push(loop);
            innerLoopsMap.set(loop, []);
        } else {
            for (let outer of outerLoops) {
                if (isPointInLoop(centroid, outer, surface, tolerance)) {
                    innerLoopsMap.get(outer).push(loop);
                    break;
                }
            }
        }
    });
    outerLoops.forEach(outerLoop => {
        const innerLoops = innerLoopsMap.get(outerLoop);
        const outerUVCurves = outerLoop.edges.map(edge => edgeToUVCurve(edge, surface, tolerance));
        for (let i = 0; i < outerUVCurves.length; i++) {
            const curr = outerUVCurves[i];
            const next = outerUVCurves[(i + 1) % outerUVCurves.length];
            const currEnd = curr.evaluate(curr.getDomain()[1]);
            next.controlPoints[0] = currEnd;
        }
        const innerUVCurves = innerLoops.map(innerLoop => {
            const curves = innerLoop.edges.map(edge => edgeToUVCurve(edge, surface, tolerance));
            for (let i = 0; i < curves.length; i++) {
                const curr = curves[i];
                const next = curves[(i + 1) % curves.length];
                const currEnd = curr.evaluate(curr.getDomain()[1]);
                next.controlPoints[0] = currEnd;
            }
            return curves;
        });
        const trimmedSurface = surface.trimWithLoop(outerUVCurves, innerUVCurves, tolerance);
        const faceLoops = [new Loop(outerLoop.edges)];
        innerLoops.forEach(innerLoop => faceLoops.push(new Loop(innerLoop.edges)));
        faces.push(new Face(trimmedSurface, faceLoops));
    });
    return faces;
}
const tests = new Tests();
const tester = new Tester(tests);
tester.runTests();