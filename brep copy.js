const defaultTolerance = 0.000001;
class Vector3D {
    constructor(x = 0, y = 0, z = 0) {
        if (typeof x !== 'number' || typeof y !== 'number' || typeof z !== 'number') {
            throw new Error('Vector3D constructor requires numeric values');
        }
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            throw new Error('Vector3D coordinates must be finite numbers');
        }
        // Store values with protection against -0
        this.x = x === 0 ? 0 : x;
        this.y = y === 0 ? 0 : y;
        this.z = z === 0 ? 0 : z;
        // Freeze the object to prevent modifications after creation
        Object.freeze(this);
    }
    add(vector) {
        if (!(vector instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        const newX = this.x + vector.x;
        const newY = this.y + vector.y;
        const newZ = this.z + vector.z;
        // Check for potential overflow/underflow
        if (!Number.isFinite(newX) || !Number.isFinite(newY) || !Number.isFinite(newZ)) {
            throw new Error('Vector addition results in non-finite values');
        }
        return new Vector3D(newX, newY, newZ);
    }
    subtract(vector) {
        if (!(vector instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        const newX = this.x - vector.x;
        const newY = this.y - vector.y;
        const newZ = this.z - vector.z;
        // Check for potential overflow/underflow
        if (!Number.isFinite(newX) || !Number.isFinite(newY) || !Number.isFinite(newZ)) {
            throw new Error('Vector subtraction results in non-finite values');
        }
        return new Vector3D(newX, newY, newZ);
    }
    multiply(scalar) {
        if (typeof scalar !== 'number') {
            throw new Error('Scalar multiplier must be a number');
        }
        if (!Number.isFinite(scalar)) {
            throw new Error('Scalar multiplier must be finite');
        }
        const newX = this.x * scalar;
        const newY = this.y * scalar;
        const newZ = this.z * scalar;
        // Check for overflow/underflow
        if (!Number.isFinite(newX) || !Number.isFinite(newY) || !Number.isFinite(newZ)) {
            throw new Error('Vector multiplication results in non-finite values');
        }
        return new Vector3D(newX, newY, newZ);
    }
    // ... other methods ...
    dot(vector) {
        if (!(vector instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        const result = this.x * vector.x + this.y * vector.y + this.z * vector.z;
        // Check for overflow/underflow
        if (!Number.isFinite(result)) {
            throw new Error('Dot product results in non-finite value');
        }
        return result;
    }
    // ... other methods ...
    cross(vector) {
        if (!(vector instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        const newX = this.y * vector.z - this.z * vector.y;
        const newY = this.z * vector.x - this.x * vector.z;
        const newZ = this.x * vector.y - this.y * vector.x;
        // Check for overflow/underflow
        if (!Number.isFinite(newX) || !Number.isFinite(newY) || !Number.isFinite(newZ)) {
            throw new Error('Cross product results in non-finite values');
        }
        return new Vector3D(newX, newY, newZ);
    }
    normalize() {
        const magnitude = Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
        // Check if vector is zero vector
        if (Math.abs(magnitude) < defaultTolerance) {
            throw new Error('Cannot normalize zero vector');
        }
        // Check for potential overflow during magnitude calculation
        if (!Number.isFinite(magnitude)) {
            throw new Error('Vector magnitude calculation results in non-finite value');
        }
        const invMag = 1 / magnitude;
        // Check for potential overflow during division
        if (!Number.isFinite(invMag)) {
            throw new Error('Vector normalization results in non-finite values');
        }
        return new Vector3D(this.x * invMag, this.y * invMag, this.z * invMag);
    }
    project(onto) {
        if (!(onto instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        // Check if projection vector is zero vector
        const ontoLengthSquared = onto.dot(onto);
        if (Math.abs(ontoLengthSquared) < defaultTolerance) {
            throw new Error('Cannot project onto zero vector');
        }
        // Calculate projection scalar
        const projectionScalar = this.dot(onto) / ontoLengthSquared;
        // Check for numerical overflow/underflow
        if (!Number.isFinite(projectionScalar)) {
            throw new Error('Projection calculation results in non-finite values');
        }
        // Return the projection vector
        return onto.multiply(projectionScalar);
    }
    reject(from) {
        if (!(from instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        // Check if rejection vector is zero vector
        const fromLengthSquared = from.dot(from);
        if (Math.abs(fromLengthSquared) < defaultTolerance) {
            throw new Error('Cannot reject from zero vector');
        }
        // Calculate projection vector
        const projection = this.project(from);
        // Rejection is the original vector minus its projection
        const rejection = this.subtract(projection);
        // Check for numerical overflow/underflow in the result
        if (!Number.isFinite(rejection.x) || !Number.isFinite(rejection.y) || !Number.isFinite(rejection.z)) {
            throw new Error('Rejection calculation results in non-finite values');
        }
        return rejection;
    }
    rotate(axis, angleRadians) {
        // Input validation
        if (!(axis instanceof Vector3D)) {
            throw new Error('Rotation axis must be a Vector3D instance');
        }
        if (typeof angleRadians !== 'number' || !Number.isFinite(angleRadians)) {
            throw new Error('Angle must be a finite number');
        }
        // Ensure axis is normalized
        const normalizedAxis = axis.normalize();
        // Calculate quaternion components
        const halfAngle = angleRadians / 2;
        const sinHalfAngle = Math.sin(halfAngle);
        // Create quaternion components [w, x, y, z]
        const qw = Math.cos(halfAngle);
        const qx = normalizedAxis.x * sinHalfAngle;
        const qy = normalizedAxis.y * sinHalfAngle;
        const qz = normalizedAxis.z * sinHalfAngle;
        // Perform quaternion rotation
        // v' = qvq*
        // where v is the vector (0, x, y, z) and q* is the conjugate quaternion
        // First multiplication (qv)
        const temp1x = qw * this.x + qy * this.z - qz * this.y;
        const temp1y = qw * this.y + qz * this.x - qx * this.z;
        const temp1z = qw * this.z + qx * this.y - qy * this.x;
        const temp1w = -qx * this.x - qy * this.y - qz * this.z;
        // Second multiplication with conjugate (qvq*)
        const newX = temp1w * -qx + temp1x * qw + temp1y * -qz - temp1z * -qy;
        const newY = temp1w * -qy + temp1y * qw + temp1z * -qx - temp1x * -qz;
        const newZ = temp1w * -qz + temp1z * qw + temp1x * -qy - temp1y * -qx;
        // Check for numerical overflow/underflow
        if (!Number.isFinite(newX) || !Number.isFinite(newY) || !Number.isFinite(newZ)) {
            throw new Error('Rotation calculation resulted in non-finite values');
        }
        // Return new rotated vector
        return new Vector3D(newX, newY, newZ);
    }
    // Add these methods to complete the implementation
    length() {
        return Math.sqrt(this.dot(this));
    }
    lengthSquared() {
        return this.dot(this);
    }
    distanceTo(vector) {
        return this.subtract(vector).length();
    }
    equals(vector, tolerance = defaultTolerance) {
        if (!(vector instanceof Vector3D)) {
            return false;
        }
        return Math.abs(this.x - vector.x) < tolerance && Math.abs(this.y - vector.y) < tolerance && Math.abs(this.z - vector.z) < tolerance;
    }
    clone() {
        return new Vector3D(this.x, this.y, this.z);
    }
}
class Matrix4x4 {
    constructor(elements) {
        // Initialize elements array with identity matrix by default
        this.elements = new Float64Array([
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            1
        ]);
        // If elements array provided, validate and set
        if (elements) {
            if (!Array.isArray(elements) && !(elements instanceof Float64Array)) {
                throw new Error('Matrix elements must be an array');
            }
            if (elements.length !== 16) {
                throw new Error('Matrix must have exactly 16 elements');
            }
            // Validate all elements are finite numbers
            if (!elements.every(e => typeof e === 'number' && Number.isFinite(e))) {
                throw new Error('All matrix elements must be finite numbers');
            }
            // Copy elements to internal array
            this.elements = new Float64Array(elements);
        }
        // Freeze the object to prevent modifications after creation
        Object.freeze(this);
    }
    multiply(matrix) {
        if (!(matrix instanceof Matrix4x4)) {
            throw new Error('Parameter must be a Matrix4x4 instance');
        }
        const result = new Float64Array(16);
        const a = this.elements;
        const b = matrix.elements;
        // Row 1
        result[0] = a[0] * b[0] + a[1] * b[4] + a[2] * b[8] + a[3] * b[12];
        result[1] = a[0] * b[1] + a[1] * b[5] + a[2] * b[9] + a[3] * b[13];
        result[2] = a[0] * b[2] + a[1] * b[6] + a[2] * b[10] + a[3] * b[14];
        result[3] = a[0] * b[3] + a[1] * b[7] + a[2] * b[11] + a[3] * b[15];
        // Row 2
        result[4] = a[4] * b[0] + a[5] * b[4] + a[6] * b[8] + a[7] * b[12];
        result[5] = a[4] * b[1] + a[5] * b[5] + a[6] * b[9] + a[7] * b[13];
        result[6] = a[4] * b[2] + a[5] * b[6] + a[6] * b[10] + a[7] * b[14];
        result[7] = a[4] * b[3] + a[5] * b[7] + a[6] * b[11] + a[7] * b[15];
        // Row 3
        result[8] = a[8] * b[0] + a[9] * b[4] + a[10] * b[8] + a[11] * b[12];
        result[9] = a[8] * b[1] + a[9] * b[5] + a[10] * b[9] + a[11] * b[13];
        result[10] = a[8] * b[2] + a[9] * b[6] + a[10] * b[10] + a[11] * b[14];
        result[11] = a[8] * b[3] + a[9] * b[7] + a[10] * b[11] + a[11] * b[15];
        // Row 4
        result[12] = a[12] * b[0] + a[13] * b[4] + a[14] * b[8] + a[15] * b[12];
        result[13] = a[12] * b[1] + a[13] * b[5] + a[14] * b[9] + a[15] * b[13];
        result[14] = a[12] * b[2] + a[13] * b[6] + a[14] * b[10] + a[15] * b[14];
        result[15] = a[12] * b[3] + a[13] * b[7] + a[14] * b[11] + a[15] * b[15];
        // Check for numerical overflow/underflow
        if (!result.every(value => Number.isFinite(value))) {
            throw new Error('Matrix multiplication resulted in non-finite values');
        }
        return new Matrix4x4(result);
    }
    inverse() {
        const m = this.elements;
        const inv = new Float64Array(16);
        // Calculate cofactors and determinant
        inv[0] = m[5] * m[10] * m[15] + m[6] * m[11] * m[13] + m[7] * m[9] * m[14] - m[5] * m[11] * m[14] - m[6] * m[9] * m[15] - m[7] * m[10] * m[13];
        inv[1] = m[1] * m[11] * m[14] + m[2] * m[9] * m[15] + m[3] * m[10] * m[13] - m[1] * m[10] * m[15] - m[2] * m[11] * m[13] - m[3] * m[9] * m[14];
        inv[2] = m[1] * m[6] * m[15] + m[2] * m[7] * m[13] + m[3] * m[5] * m[14] - m[1] * m[7] * m[14] - m[2] * m[5] * m[15] - m[3] * m[6] * m[13];
        inv[3] = m[1] * m[7] * m[10] + m[2] * m[5] * m[11] + m[3] * m[6] * m[9] - m[1] * m[6] * m[11] - m[2] * m[7] * m[9] - m[3] * m[5] * m[10];
        inv[4] = m[4] * m[11] * m[14] + m[6] * m[8] * m[15] + m[7] * m[10] * m[12] - m[4] * m[10] * m[15] - m[6] * m[11] * m[12] - m[7] * m[8] * m[14];
        inv[5] = m[0] * m[10] * m[15] + m[2] * m[11] * m[12] + m[3] * m[8] * m[14] - m[0] * m[11] * m[14] - m[2] * m[8] * m[15] - m[3] * m[10] * m[12];
        inv[6] = m[0] * m[7] * m[14] + m[2] * m[4] * m[15] + m[3] * m[6] * m[12] - m[0] * m[6] * m[15] - m[2] * m[7] * m[12] - m[3] * m[4] * m[14];
        inv[7] = m[0] * m[6] * m[11] + m[2] * m[7] * m[8] + m[3] * m[4] * m[10] - m[0] * m[7] * m[10] - m[2] * m[4] * m[11] - m[3] * m[6] * m[8];
        inv[8] = m[4] * m[9] * m[15] + m[5] * m[11] * m[12] + m[7] * m[8] * m[13] - m[4] * m[11] * m[13] - m[5] * m[8] * m[15] - m[7] * m[9] * m[12];
        inv[9] = m[0] * m[11] * m[13] + m[1] * m[8] * m[15] + m[3] * m[9] * m[12] - m[0] * m[9] * m[15] - m[1] * m[11] * m[12] - m[3] * m[8] * m[13];
        inv[10] = m[0] * m[5] * m[15] + m[1] * m[7] * m[12] + m[3] * m[4] * m[13] - m[0] * m[7] * m[13] - m[1] * m[4] * m[15] - m[3] * m[5] * m[12];
        inv[11] = m[0] * m[7] * m[9] + m[1] * m[4] * m[11] + m[3] * m[5] * m[8] - m[0] * m[5] * m[11] - m[1] * m[7] * m[8] - m[3] * m[4] * m[9];
        inv[12] = m[4] * m[10] * m[13] + m[5] * m[8] * m[14] + m[6] * m[9] * m[12] - m[4] * m[9] * m[14] - m[5] * m[10] * m[12] - m[6] * m[8] * m[13];
        inv[13] = m[0] * m[9] * m[14] + m[1] * m[10] * m[12] + m[2] * m[8] * m[13] - m[0] * m[10] * m[13] - m[1] * m[8] * m[14] - m[2] * m[9] * m[12];
        inv[14] = m[0] * m[6] * m[13] + m[1] * m[4] * m[14] + m[2] * m[5] * m[12] - m[0] * m[5] * m[14] - m[1] * m[6] * m[12] - m[2] * m[4] * m[13];
        inv[15] = m[0] * m[5] * m[10] + m[1] * m[6] * m[8] + m[2] * m[4] * m[9] - m[0] * m[6] * m[9] - m[1] * m[4] * m[10] - m[2] * m[5] * m[8];
        // Calculate determinant using first row of cofactors
        const det = m[0] * inv[0] - m[1] * inv[4] + m[2] * inv[8] - m[3] * inv[12];
        // Check if matrix is invertible
        if (Math.abs(det) < Number.EPSILON) {
            throw new Error('Matrix is not invertible');
        }
        // Multiply adjugate matrix by 1/determinant
        const invDet = 1 / det;
        for (let i = 0; i < 16; i++) {
            inv[i] *= invDet;
            // Check for numerical overflow/underflow
            if (!Number.isFinite(inv[i])) {
                throw new Error('Matrix inversion resulted in non-finite values');
            }
        }
        return new Matrix4x4(inv);
    }
    transpose() {
        const m = this.elements;
        const result = new Float64Array([
            m[0],
            m[4],
            m[8],
            m[12],
            m[1],
            m[5],
            m[9],
            m[13],
            m[2],
            m[6],
            m[10],
            m[14],
            m[3],
            m[7],
            m[11],
            m[15]
        ]);
        // Check for numerical validity
        if (!result.every(value => Number.isFinite(value))) {
            throw new Error('Matrix transposition resulted in non-finite values');
        }
        return new Matrix4x4(result);
    }
    createTransformation(translation = {
        x: 0,
        y: 0,
        z: 0
    }, rotation = {
        x: 0,
        y: 0,
        z: 0
    }, scale = {
        x: 1,
        y: 1,
        z: 1
    }) {
        // Validate input parameters
        const params = [
            translation.x,
            translation.y,
            translation.z,
            rotation.x,
            rotation.y,
            rotation.z,
            scale.x,
            scale.y,
            scale.z
        ];
        if (!params.every(param => typeof param === 'number' && Number.isFinite(param))) {
            throw new Error('All transformation parameters must be finite numbers');
        }
        // Create rotation matrices for each axis
        const cosX = Math.cos(rotation.x);
        const sinX = Math.sin(rotation.x);
        const cosY = Math.cos(rotation.y);
        const sinY = Math.sin(rotation.y);
        const cosZ = Math.cos(rotation.z);
        const sinZ = Math.sin(rotation.z);
        // Create rotation matrix around X axis
        const rotX = new Matrix4x4([
            1,
            0,
            0,
            0,
            0,
            cosX,
            -sinX,
            0,
            0,
            sinX,
            cosX,
            0,
            0,
            0,
            0,
            1
        ]);
        // Create rotation matrix around Y axis
        const rotY = new Matrix4x4([
            cosY,
            0,
            sinY,
            0,
            0,
            1,
            0,
            0,
            -sinY,
            0,
            cosY,
            0,
            0,
            0,
            0,
            1
        ]);
        // Create rotation matrix around Z axis
        const rotZ = new Matrix4x4([
            cosZ,
            -sinZ,
            0,
            0,
            sinZ,
            cosZ,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            1
        ]);
        // Create scale matrix
        const scaleMatrix = new Matrix4x4([
            scale.x,
            0,
            0,
            0,
            0,
            scale.y,
            0,
            0,
            0,
            0,
            scale.z,
            0,
            0,
            0,
            0,
            1
        ]);
        // Create translation matrix
        const translationMatrix = new Matrix4x4([
            1,
            0,
            0,
            translation.x,
            0,
            1,
            0,
            translation.y,
            0,
            0,
            1,
            translation.z,
            0,
            0,
            0,
            1
        ]);
        try {
            // Combine transformations: Translation * RotZ * RotY * RotX * Scale
            const rotationMatrix = rotZ.multiply(rotY).multiply(rotX);
            const transformationMatrix = translationMatrix.multiply(rotationMatrix).multiply(scaleMatrix);
            // Validate the final result
            if (!transformationMatrix.elements.every(value => Number.isFinite(value))) {
                throw new Error('Transformation resulted in non-finite values');
            }
            return transformationMatrix;
        } catch (error) {
            throw new Error(`Failed to create transformation matrix: ${error.message}`);
        }
    }
    // ... other methods ...
    determinant() {
        const m = this.elements;
        // Calculate cofactors for the first row
        const c00 = m[5] * (m[10] * m[15] - m[11] * m[14]) - m[9] * (m[6] * m[15] - m[7] * m[14]) + m[13] * (m[6] * m[11] - m[7] * m[10]);
        const c01 = m[4] * (m[10] * m[15] - m[11] * m[14]) - m[8] * (m[6] * m[15] - m[7] * m[14]) + m[12] * (m[6] * m[11] - m[7] * m[10]);
        const c02 = m[4] * (m[9] * m[15] - m[11] * m[13]) - m[8] * (m[5] * m[15] - m[7] * m[13]) + m[12] * (m[5] * m[11] - m[7] * m[9]);
        const c03 = m[4] * (m[9] * m[14] - m[10] * m[13]) - m[8] * (m[5] * m[14] - m[6] * m[13]) + m[12] * (m[5] * m[10] - m[6] * m[9]);
        // Calculate determinant using the first row and its cofactors
        const det = m[0] * c00 - m[1] * c01 + m[2] * c02 - m[3] * c03;
        // Check for numerical validity
        if (!Number.isFinite(det)) {
            throw new Error('Determinant calculation resulted in non-finite value');
        }
        // Handle -0 case
        return det === 0 ? 0 : det;
    }
}
class KnotVector {
    constructor(knots = []) {
        // Input validation
        if (!Array.isArray(knots)) {
            throw new Error('Knot vector must be an array');
        }
        // Validate each knot value is a finite number
        if (!knots.every(knot => typeof knot === 'number' && Number.isFinite(knot))) {
            throw new Error('All knots must be finite numbers');
        }
        // Create a copy of the knots array using Float64Array for numerical precision
        this.knots = new Float64Array(knots);
        // Validate monotonic non-decreasing sequence
        for (let i = 1; i < this.knots.length; i++) {
            if (this.knots[i] < this.knots[i - 1]) {
                throw new Error('Knot vector must be monotonically non-decreasing');
            }
        }
        // Store additional properties
        this.length = this.knots.length;
        this.domain = {
            min: this.knots[0],
            max: this.knots[this.length - 1]
        };
        // Calculate multiplicity array
        this.multiplicities = this.calculateMultiplicities();
        // Freeze the object to prevent modifications after creation
        Object.freeze(this.knots);
        Object.freeze(this.domain);
        Object.freeze(this.multiplicities);
        Object.freeze(this);
    }
    normalize() {
        // If knot vector is empty or has only one knot, return a new empty or single-knot vector
        if (this.knots.length <= 1) {
            return new KnotVector([...this.knots]);
        }
        // Get the domain range
        const min = this.domain.min;
        const max = this.domain.max;
        const range = max - min;
        // Check if the knot vector is already normalized
        if (Math.abs(min) < Number.EPSILON && Math.abs(max - 1) < Number.EPSILON) {
            return new KnotVector([...this.knots]);
        }
        // Check for zero range (all knots are the same value)
        if (Math.abs(range) < Number.EPSILON) {
            throw new Error('Cannot normalize knot vector with zero range');
        }
        // Create new normalized knot vector
        const normalizedKnots = new Float64Array(this.knots.length);
        // Normalize each knot value to [0,1] interval
        for (let i = 0; i < this.knots.length; i++) {
            const normalizedValue = (this.knots[i] - min) / range;
            // Handle potential numerical precision issues
            if (normalizedValue < 0 || normalizedValue > 1) {
                throw new Error('Normalization resulted in values outside [0,1] range');
            }
            // Ensure exact 0 and 1 for the endpoints
            if (Math.abs(normalizedValue) < Number.EPSILON) {
                normalizedKnots[i] = 0;
            } else if (Math.abs(normalizedValue - 1) < Number.EPSILON) {
                normalizedKnots[i] = 1;
            } else {
                normalizedKnots[i] = normalizedValue;
            }
        }
        // Verify the normalization maintains monotonicity
        for (let i = 1; i < normalizedKnots.length; i++) {
            if (normalizedKnots[i] < normalizedKnots[i - 1]) {
                throw new Error('Normalization failed to maintain monotonicity');
            }
        }
        // Return new KnotVector instance with normalized knots
        return new KnotVector(normalizedKnots);
    }
    insert(u, multiplicity = 1) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Knot value must be a finite number');
        }
        if (!Number.isInteger(multiplicity) || multiplicity < 1) {
            throw new Error('Multiplicity must be a positive integer');
        }
        // Validate knot value is within domain
        if (u < this.domain.min || u > this.domain.max) {
            throw new Error('Knot value must be within the knot vector domain');
        }
        // Find insertion position
        let insertionIndex = this.knots.findIndex(knot => knot > u);
        if (insertionIndex === -1) {
            insertionIndex = this.knots.length;
        }
        // Calculate current multiplicity at u
        const currentMultiplicity = this.multiplicities.get(u) || 0;
        const actualMultiplicity = Math.min(multiplicity, this.degree + 1 - currentMultiplicity);
        if (actualMultiplicity <= 0) {
            return new KnotVector([...this.knots]);
        }
        // Create new knot array with inserted knots
        const newKnots = new Float64Array(this.knots.length + actualMultiplicity);
        // Copy knots before insertion point
        newKnots.set(this.knots.slice(0, insertionIndex));
        // Insert new knots
        for (let i = 0; i < actualMultiplicity; i++) {
            newKnots[insertionIndex + i] = u;
        }
        // Copy remaining knots
        newKnots.set(this.knots.slice(insertionIndex), insertionIndex + actualMultiplicity);
        // Validate new knot vector maintains monotonicity
        for (let i = 1; i < newKnots.length; i++) {
            if (newKnots[i] < newKnots[i - 1]) {
                throw new Error('Knot insertion resulted in non-monotonic sequence');
            }
        }
        try {
            // Return new KnotVector instance with inserted knots
            return new KnotVector(newKnots);
        } catch (error) {
            throw new Error(`Failed to create new knot vector: ${error.message}`);
        }
    }
    remove(u, multiplicity = 1) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Knot value must be a finite number');
        }
        if (!Number.isInteger(multiplicity) || multiplicity < 1) {
            throw new Error('Multiplicity must be a positive integer');
        }
        // Check if knot exists in vector
        const currentMultiplicity = this.multiplicities.get(u);
        if (!currentMultiplicity) {
            return new KnotVector([...this.knots]);
        }
        // Return copy if knot not found
        // Calculate how many knots we can actually remove
        const actualMultiplicity = Math.min(multiplicity, currentMultiplicity - this.degree + 1);
        if (actualMultiplicity <= 0) {
            return new KnotVector([...this.knots]);
        }
        // Return copy if we can't remove any
        // Find the start index of the knot to remove
        let startIndex = -1;
        for (let i = 0; i < this.knots.length; i++) {
            if (Math.abs(this.knots[i] - u) < Number.EPSILON) {
                startIndex = i;
                break;
            }
        }
        if (startIndex === -1) {
            throw new Error('Internal error: Knot not found despite having multiplicity');
        }
        // Create new knot array without removed knots
        const newKnots = new Float64Array(this.knots.length - actualMultiplicity);
        // Copy knots before removal point
        newKnots.set(this.knots.slice(0, startIndex));
        // Copy remaining knots after skipping removed ones
        const remainingKnots = this.knots.slice(startIndex + actualMultiplicity);
        newKnots.set(remainingKnots, startIndex);
        // Validate new knot vector
        if (newKnots.length < 2) {
            throw new Error('Cannot remove knots: Would result in invalid knot vector');
        }
        // Validate monotonicity is maintained
        for (let i = 1; i < newKnots.length; i++) {
            if (newKnots[i] < newKnots[i - 1]) {
                throw new Error('Knot removal would violate monotonicity');
            }
        }
        try {
            // Return new KnotVector instance with removed knots
            return new KnotVector(newKnots);
        } catch (error) {
            throw new Error(`Failed to create new knot vector: ${error.message}`);
        }
    }
    refine(refinementLevel = 1) {
        // Input validation
        if (!Number.isInteger(refinementLevel) || refinementLevel < 1) {
            throw new Error('Refinement level must be a positive integer');
        }
        // If knot vector is empty or has only one knot, return a copy
        if (this.knots.length <= 1) {
            return new KnotVector([...this.knots]);
        }
        // Calculate new knots to be inserted
        const newKnots = new Set();
        // For each refinement level
        for (let level = 0; level < refinementLevel; level++) {
            for (let i = 0; i < this.knots.length - 1; i++) {
                const currentKnot = this.knots[i];
                const nextKnot = this.knots[i + 1];
                // Only insert new knot if the knots are different
                if (Math.abs(nextKnot - currentKnot) > Number.EPSILON) {
                    // Calculate midpoint
                    const midKnot = (currentKnot + nextKnot) / 2;
                    // Check for numerical validity
                    if (!Number.isFinite(midKnot)) {
                        throw new Error('Refinement resulted in non-finite knot value');
                    }
                    // Add to new knots set
                    newKnots.add(midKnot);
                }
            }
        }
        // Convert set to array and sort
        const knotsToInsert = Array.from(newKnots).sort((a, b) => a - b);
        // Create new array with sufficient size
        const refinedKnots = new Float64Array(this.knots.length + knotsToInsert.length);
        // Merge original and new knots while maintaining order
        let originalIndex = 0;
        let newIndex = 0;
        let resultIndex = 0;
        while (originalIndex < this.knots.length && newIndex < knotsToInsert.length) {
            if (this.knots[originalIndex] <= knotsToInsert[newIndex]) {
                refinedKnots[resultIndex] = this.knots[originalIndex];
                originalIndex++;
            } else {
                refinedKnots[resultIndex] = knotsToInsert[newIndex];
                newIndex++;
            }
            resultIndex++;
        }
        // Copy any remaining original knots
        while (originalIndex < this.knots.length) {
            refinedKnots[resultIndex] = this.knots[originalIndex];
            originalIndex++;
            resultIndex++;
        }
        // Copy any remaining new knots
        while (newIndex < knotsToInsert.length) {
            refinedKnots[resultIndex] = knotsToInsert[newIndex];
            newIndex++;
            resultIndex++;
        }
        // Validate monotonicity of refined knot vector
        for (let i = 1; i < refinedKnots.length; i++) {
            if (refinedKnots[i] < refinedKnots[i - 1]) {
                throw new Error('Refinement resulted in non-monotonic knot vector');
            }
        }
        try {
            // Create and return new KnotVector instance with refined knots
            return new KnotVector(refinedKnots);
        } catch (error) {
            throw new Error(`Failed to create refined knot vector: ${error.message}`);
        }
    }
    validate() {
        // Check if knot vector exists and has valid length
        if (!this.knots || this.knots.length < 2) {
            throw new Error('Knot vector must contain at least 2 values');
        }
        // Check data type and finiteness of each knot
        for (let i = 0; i < this.knots.length; i++) {
            const knot = this.knots[i];
            if (typeof knot !== 'number') {
                throw new Error(`Invalid knot type at index ${i}: must be a number`);
            }
            if (!Number.isFinite(knot)) {
                throw new Error(`Invalid knot value at index ${i}: must be finite`);
            }
        }
        // Check monotonicity (non-decreasing sequence)
        for (let i = 1; i < this.knots.length; i++) {
            if (this.knots[i] < this.knots[i - 1]) {
                throw new Error(`Knot vector violates monotonicity at index ${i}`);
            }
        }
        // Check for valid multiplicity (cannot exceed degree + 1)
        if (this.degree !== undefined) {
            const maxMultiplicity = this.degree + 1;
            let currentKnot = this.knots[0];
            let currentMultiplicity = 1;
            for (let i = 1; i < this.knots.length; i++) {
                if (Math.abs(this.knots[i] - currentKnot) < Number.EPSILON) {
                    currentMultiplicity++;
                    if (currentMultiplicity > maxMultiplicity) {
                        throw new Error(`Knot multiplicity exceeds degree + 1 at value ${currentKnot}`);
                    }
                } else {
                    currentKnot = this.knots[i];
                    currentMultiplicity = 1;
                }
            }
        }
        // Check domain validity
        if (this.domain) {
            if (!Number.isFinite(this.domain.min) || !Number.isFinite(this.domain.max)) {
                throw new Error('Invalid domain bounds');
            }
            if (this.domain.max < this.domain.min) {
                throw new Error('Invalid domain range: max must be greater than or equal to min');
            }
            if (Math.abs(this.domain.min - this.knots[0]) > Number.EPSILON || Math.abs(this.domain.max - this.knots[this.knots.length - 1]) > Number.EPSILON) {
                throw new Error('Domain bounds do not match knot vector endpoints');
            }
        }
        // Check multiplicity map validity if it exists
        if (this.multiplicities) {
            for (const [knot, multiplicity] of this.multiplicities) {
                if (!Number.isFinite(knot) || !Number.isInteger(multiplicity) || multiplicity < 1) {
                    throw new Error('Invalid multiplicity map entry');
                }
            }
        }
        // Return true if all validations pass
        return true;
    }
    calculateMultiplicities() {
        const multiplicities = new Map();
        let currentKnot = this.knots[0];
        let count = 1;
        for (let i = 1; i < this.knots.length; i++) {
            if (Math.abs(this.knots[i] - currentKnot) < Number.EPSILON) {
                count++;
            } else {
                multiplicities.set(currentKnot, count);
                currentKnot = this.knots[i];
                count = 1;
            }
        }
        multiplicities.set(currentKnot, count);
        return multiplicities;
    }
}
class ControlPoint {
    constructor(x = 0, y = 0, z = 0, w = 1) {
        // Validate input types
        if (typeof x !== 'number' || typeof y !== 'number' || typeof z !== 'number' || typeof w !== 'number') {
            throw new Error('Control point coordinates must be numbers');
        }
        // Validate for finite values
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z) || !Number.isFinite(w)) {
            throw new Error('Control point coordinates must be finite numbers');
        }
        // Validate weight (w) is non-zero
        if (Math.abs(w) < Number.EPSILON) {
            throw new Error('Weight (w) must be non-zero');
        }
        // Store coordinates, protecting against -0
        this.x = x === 0 ? 0 : x;
        this.y = y === 0 ? 0 : y;
        this.z = z === 0 ? 0 : z;
        this.w = w === 0 ? 1 : w;
        // Freeze the object to prevent modifications after creation
        Object.freeze(this);
    }
    toHomogeneous() {
        // Since the control point is already stored in homogeneous coordinates,
        // we return a new array with the coordinates
        if (!Number.isFinite(this.x) || !Number.isFinite(this.y) || !Number.isFinite(this.z) || !Number.isFinite(this.w)) {
            throw new Error('Invalid control point coordinates');
        }
        // Return array of homogeneous coordinates [x, y, z, w]
        return Object.freeze([
            this.x === 0 ? 0 : this.x,
            // Protect against -0
            this.y === 0 ? 0 : this.y,
            this.z === 0 ? 0 : this.z,
            this.w === 0 ? 1 : this.w
        ]);
    }
    fromHomogeneous(coords) {
        // Input validation
        if (!Array.isArray(coords) && !(coords instanceof Float64Array)) {
            throw new Error('Homogeneous coordinates must be provided as an array');
        }
        if (coords.length !== 4) {
            throw new Error('Homogeneous coordinates must have exactly 4 components');
        }
        // Validate each coordinate is a finite number
        if (!coords.every(coord => typeof coord === 'number' && Number.isFinite(coord))) {
            throw new Error('All coordinates must be finite numbers');
        }
        const [x, y, z, w] = coords;
        // Validate w is non-zero
        if (Math.abs(w) < Number.EPSILON) {
            throw new Error('W coordinate must be non-zero');
        }
        try {
            // Create and return a new ControlPoint instance
            return new ControlPoint(x, y, z, w);
        } catch (error) {
            throw new Error(`Failed to create control point: ${error.message}`);
        }
    }
    weight() {
        // Since w is stored and immutable, we can simply return it
        // No need for validation since constructor ensures w is valid
        return this.w;
    }
    // ... other methods ...
    position() {
        // Check if the control point coordinates are valid
        if (!Number.isFinite(this.x) || !Number.isFinite(this.y) || !Number.isFinite(this.z) || !Number.isFinite(this.w)) {
            throw new Error('Invalid control point coordinates');
        }
        // Check for zero weight
        if (Math.abs(this.w) < Number.EPSILON) {
            throw new Error('Cannot compute position: weight is zero');
        }
        // Perform perspective division to convert from homogeneous to Cartesian coordinates
        const invW = 1 / this.w;
        const x = this.x * invW;
        const y = this.y * invW;
        const z = this.z * invW;
        // Check for numerical overflow/underflow in results
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            throw new Error('Position calculation resulted in non-finite values');
        }
        // Return frozen array of position coordinates [x, y, z]
        return Object.freeze([
            x === 0 ? 0 : x,
            // Protect against -0
            y === 0 ? 0 : y,
            z === 0 ? 0 : z
        ]);
    }
}
class BasisFunction {
    constructor() {
        this.degree = 0;
        this.knotVector = null;
    }
    evaluate(i, p, u) {
        // Input validation
        if (!Number.isInteger(i) || i < 0) {
            throw new Error('Index i must be a non-negative integer');
        }
        if (!Number.isInteger(p) || p < 0) {
            throw new Error('Degree p must be a non-negative integer');
        }
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (!this.knotVector || !(this.knotVector instanceof KnotVector)) {
            throw new Error('Valid knot vector not initialized');
        }
        if (i >= this.knotVector.length - p) {
            throw new Error('Index i is out of valid range for given degree');
        }
        // Get knot vector values
        const knots = this.knotVector.knots;
        // Check if u is within the domain
        if (u < knots[0] || u > knots[knots.length - 1]) {
            throw new Error('Parameter u is outside the knot vector domain');
        }
        // Special case: zero degree basis function
        if (p === 0) {
            return u >= knots[i] && u < knots[i + 1] ? 1 : 0;
        }
        // Initialize storage for basis function values
        const N = new Array(p + 1).fill(0);
        // Calculate zeroth degree basis functions
        for (let j = 0; j <= p; j++) {
            N[j] = u >= knots[i + j] && u < knots[i + j + 1] ? 1 : 0;
        }
        // Apply Cox-de Boor recursion formula
        for (let k = 1; k <= p; k++) {
            for (let j = 0; j < p - k + 1; j++) {
                // Calculate left term
                let leftTerm = 0;
                const leftDenom = knots[i + j + k] - knots[i + j];
                if (Math.abs(leftDenom) > Number.EPSILON) {
                    leftTerm = (u - knots[i + j]) * N[j] / leftDenom;
                }
                // Calculate right term
                let rightTerm = 0;
                const rightDenom = knots[i + j + k + 1] - knots[i + j + 1];
                if (Math.abs(rightDenom) > Number.EPSILON) {
                    rightTerm = (knots[i + j + k + 1] - u) * N[j + 1] / rightDenom;
                }
                // Combine terms
                N[j] = leftTerm + rightTerm;
                // Check for numerical validity
                if (!Number.isFinite(N[j])) {
                    throw new Error('Basis function evaluation resulted in non-finite value');
                }
            }
        }
        // Return the final basis function value
        return N[0];
    }
    evaluateDerivative(i, p, u, derivative = 1) {
        // Input validation
        if (!Number.isInteger(i) || i < 0) {
            throw new Error('Index i must be a non-negative integer');
        }
        if (!Number.isInteger(p) || p < 0) {
            throw new Error('Degree p must be a non-negative integer');
        }
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (!Number.isInteger(derivative) || derivative < 0) {
            throw new Error('Derivative order must be a non-negative integer');
        }
        if (!this.knotVector || !(this.knotVector instanceof KnotVector)) {
            throw new Error('Valid knot vector not initialized');
        }
        if (i >= this.knotVector.length - p) {
            throw new Error('Index i is out of valid range for given degree');
        }
        if (derivative > p) {
            return 0;
        }
        // All derivatives higher than p are zero
        const knots = this.knotVector.knots;
        // Check if u is within the domain
        if (u < knots[0] || u > knots[knots.length - 1]) {
            throw new Error('Parameter u is outside the knot vector domain');
        }
        // Special case: zero derivative
        if (derivative === 0) {
            return this.evaluate(i, p, u);
        }
        try {
            // Calculate the derivative using the derivative formula
            let result = p / (knots[i + p] - knots[i]) * this.evaluate(i, p - 1, u);
            if (i + 1 <= knots.length - p) {
                result -= p / (knots[i + p + 1] - knots[i + 1]) * this.evaluate(i + 1, p - 1, u);
            }
            // For higher order derivatives, recursively calculate
            if (derivative > 1) {
                const factor = p / (knots[i + p] - knots[i]);
                if (Number.isFinite(factor)) {
                    result = factor * this.evaluateDerivative(i, p - 1, u, derivative - 1);
                }
                if (i + 1 <= knots.length - p) {
                    const nextFactor = p / (knots[i + p + 1] - knots[i + 1]);
                    if (Number.isFinite(nextFactor)) {
                        result -= nextFactor * this.evaluateDerivative(i + 1, p - 1, u, derivative - 1);
                    }
                }
            }
            // Check for numerical validity
            if (!Number.isFinite(result)) {
                throw new Error('Derivative calculation resulted in non-finite value');
            }
            return result;
        } catch (error) {
            throw new Error(`Failed to calculate derivative: ${error.message}`);
        }
    }
    evaluateAll(i, p, u, derivativeOrder = 0) {
        // Input validation
        if (!Number.isInteger(i) || i < 0) {
            throw new Error('Index i must be a non-negative integer');
        }
        if (!Number.isInteger(p) || p < 0) {
            throw new Error('Degree p must be a non-negative integer');
        }
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (!Number.isInteger(derivativeOrder) || derivativeOrder < 0) {
            throw new Error('Derivative order must be a non-negative integer');
        }
        if (!this.knotVector || !(this.knotVector instanceof KnotVector)) {
            throw new Error('Valid knot vector not initialized');
        }
        if (i >= this.knotVector.length - p) {
            throw new Error('Index i is out of valid range for given degree');
        }
        const knots = this.knotVector.knots;
        const n = p + 1;
        // Number of basis functions
        const ndu = new Array(n).fill(0).map(() => new Array(n).fill(0));
        const left = new Array(n).fill(0);
        const right = new Array(n).fill(0);
        const derivatives = new Array(derivativeOrder + 1).fill(0).map(() => new Array(n).fill(0));
        // Initialize zeroth degree basis function
        ndu[0][0] = 1;
        // Compute triangular table of basis functions
        for (let j = 1; j < n; j++) {
            left[j] = u - knots[i + 1 - j];
            right[j] = knots[i + j] - u;
            let saved = 0;
            for (let r = 0; r < j; r++) {
                // Lower triangle
                ndu[j][r] = right[r + 1] + left[j - r];
                const temp = ndu[r][j - 1] / ndu[j][r];
                // Upper triangle
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }
        // Load the basis functions
        for (let j = 0; j <= p; j++) {
            derivatives[0][j] = ndu[j][p];
        }
        // Calculate derivatives if requested
        if (derivativeOrder > 0) {
            // Compute derivatives using triangular table
            for (let r = 0; r <= p; r++) {
                let s1 = 0;
                let s2 = 1;
                const a = new Array(2).fill(0).map(() => new Array(n).fill(0));
                // Loop for derivative orders
                for (let k = 1; k <= derivativeOrder; k++) {
                    let d = 0;
                    const rk = r - k;
                    const pk = p - k;
                    if (r >= k) {
                        a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                        d = a[s2][0] * ndu[rk][pk];
                    }
                    const j1 = rk >= -1 ? 1 : -rk;
                    const j2 = r - 1 <= pk ? k - 1 : p - r;
                    for (let j = j1; j <= j2; j++) {
                        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                        d += a[s2][j] * ndu[rk + j][pk];
                    }
                    if (r <= pk) {
                        a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                        d += a[s2][k] * ndu[r][pk];
                    }
                    derivatives[k][r] = d;
                    // Switch indices for next iteration
                    const temp = s1;
                    s1 = s2;
                    s2 = temp;
                }
            }
            // Multiply through by correct factors
            let factor = p;
            for (let k = 1; k <= derivativeOrder; k++) {
                for (let j = 0; j <= p; j++) {
                    derivatives[k][j] *= factor;
                }
                factor *= p - k;
            }
        }
        // Check for numerical validity
        for (const row of derivatives) {
            for (const value of row) {
                if (!Number.isFinite(value)) {
                    throw new Error('Basis function evaluation resulted in non-finite values');
                }
            }
        }
        return derivatives;
    }
    supportInterval(i, p) {
        // Input validation
        if (!Number.isInteger(i) || i < 0) {
            throw new Error('Index i must be a non-negative integer');
        }
        if (!Number.isInteger(p) || p < 0) {
            throw new Error('Degree p must be a non-negative integer');
        }
        if (!this.knotVector || !(this.knotVector instanceof KnotVector)) {
            throw new Error('Valid knot vector not initialized');
        }
        if (i >= this.knotVector.length - p) {
            throw new Error('Index i is out of valid range for given degree');
        }
        const knots = this.knotVector.knots;
        // The support interval for basis function N_{i,p} is [u_i, u_{i+p+1})
        const start = knots[i];
        const end = knots[i + p + 1];
        // Validate interval bounds
        if (!Number.isFinite(start) || !Number.isFinite(end)) {
            throw new Error('Support interval bounds are not finite');
        }
        if (end < start) {
            throw new Error('Invalid support interval: end is less than start');
        }
        // Return interval as an object with start and end properties
        return Object.freeze({
            start: start,
            end: end,
            isEmpty: Math.abs(end - start) < Number.EPSILON,
            contains: function (u) {
                return u >= this.start && u < this.end;
            }
        });
    }
}
class NurbsCurve {
    constructor(controlPoints = [], knotVector = [], degree = 3) {
        // Validate degree
        if (!Number.isInteger(degree) || degree < 1) {
            throw new Error('Degree must be a positive integer');
        }
        // Validate control points
        if (!Array.isArray(controlPoints)) {
            throw new Error('Control points must be provided as an array');
        }
        // Verify all control points are valid ControlPoint instances
        if (!controlPoints.every(point => point instanceof ControlPoint)) {
            throw new Error('All control points must be instances of ControlPoint class');
        }
        // Validate knot vector
        if (!(knotVector instanceof KnotVector)) {
            if (Array.isArray(knotVector)) {
                try {
                    knotVector = new KnotVector(knotVector);
                } catch (error) {
                    throw new Error(`Invalid knot vector: ${error.message}`);
                }
            } else {
                throw new Error('Knot vector must be an instance of KnotVector or an array');
            }
        }
        // Validate relationship between control points, degree, and knot vector
        const expectedKnotLength = controlPoints.length + degree + 1;
        if (knotVector.length !== expectedKnotLength) {
            throw new Error(`Knot vector length must be ${expectedKnotLength} for ${controlPoints.length} control points of degree ${degree}`);
        }
        // Minimum number of control points must be degree + 1
        if (controlPoints.length < degree + 1) {
            throw new Error(`At least ${degree + 1} control points are required for degree ${degree}`);
        }
        // Store properties using private references to prevent modification
        this.controlPoints = Object.freeze([...controlPoints]);
        this.knotVector = knotVector;
        this.degree = degree;
        // Calculate and store curve domain
        this.domain = Object.freeze({
            min: knotVector.knots[degree],
            max: knotVector.knots[knotVector.length - degree - 1]
        });
        // Initialize basis function calculator
        this.basisFunction = new BasisFunction();
        this.basisFunction.knotVector = this.knotVector;
        this.basisFunction.degree = this.degree;
        // Freeze the entire object to prevent modifications after creation
        Object.freeze(this);
    }
    // ... other methods ...
    evaluate(u) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        // Check if u is within the domain
        if (u < this.domain.min || u > this.domain.max) {
            throw new Error('Parameter u is outside the curve domain');
        }
        // Initialize numerator and denominator for rational calculation
        let numeratorX = 0;
        let numeratorY = 0;
        let numeratorZ = 0;
        let denominator = 0;
        // Calculate the sum using basis functions and control points
        for (let i = 0; i < this.controlPoints.length; i++) {
            // Calculate basis function value
            const basisValue = this.basisFunction.evaluate(i, this.degree, u);
            // Get control point weight and coordinates
            const controlPoint = this.controlPoints[i];
            const weight = controlPoint.weight();
            const [x, y, z] = controlPoint.position();
            // Calculate weighted contribution
            const weightedBasis = basisValue * weight;
            // Add to numerator and denominator
            numeratorX += weightedBasis * x;
            numeratorY += weightedBasis * y;
            numeratorZ += weightedBasis * z;
            denominator += weightedBasis;
        }
        // Check for zero denominator
        if (Math.abs(denominator) < Number.EPSILON) {
            throw new Error('Invalid curve evaluation: zero denominator');
        }
        // Calculate final point coordinates
        const invDenom = 1 / denominator;
        const x = numeratorX * invDenom;
        const y = numeratorY * invDenom;
        const z = numeratorZ * invDenom;
        // Check for numerical validity
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            throw new Error('Curve evaluation resulted in non-finite values');
        }
        // Return point on curve as Vector3D instance
        return new Vector3D(x === 0 ? 0 : x, // Protect against -0
            y === 0 ? 0 : y, z === 0 ? 0 : z);
    }
    derivative(u, order = 1) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (!Number.isInteger(order) || order < 1) {
            throw new Error('Derivative order must be a positive integer');
        }
        if (order > this.degree) {
            // Higher order derivatives than the degree are zero vectors
            return new Vector3D(0, 0, 0);
        }
        // Check if u is within the domain
        if (u < this.domain.min || u > this.domain.max) {
            throw new Error('Parameter u is outside the curve domain');
        }
        // Calculate basis functions and their derivatives up to required order
        const basisDerivatives = [];
        for (let i = 0; i < this.controlPoints.length; i++) {
            basisDerivatives.push(this.basisFunction.evaluateAll(i, this.degree, u, order));
        }
        // Initialize arrays for storing weighted sums
        const A = new Array(order + 1).fill(0).map(() => ({
            x: 0,
            y: 0,
            z: 0
        }));
        let w = 0;
        const wDerivatives = new Array(order).fill(0);
        // Calculate weighted sums for each derivative order
        for (let i = 0; i < this.controlPoints.length; i++) {
            const controlPoint = this.controlPoints[i];
            const weight = controlPoint.weight();
            const [px, py, pz] = controlPoint.position();
            // Calculate weighted contributions for each order
            for (let k = 0; k <= order; k++) {
                const basis = basisDerivatives[i][k];
                const weightedBasis = basis * weight;
                if (k === 0) {
                    w += weightedBasis;
                    A[0].x += weightedBasis * px;
                    A[0].y += weightedBasis * py;
                    A[0].z += weightedBasis * pz;
                } else {
                    wDerivatives[k - 1] += weightedBasis;
                    A[k].x += weightedBasis * px;
                    A[k].y += weightedBasis * py;
                    A[k].z += weightedBasis * pz;
                }
            }
        }
        // Apply quotient rule to compute the actual derivative
        const result = {
            x: 0,
            y: 0,
            z: 0
        };
        // Using generalized Leibniz rule for quotient rule
        for (let k = 0; k <= order; k++) {
            const sign = k % 2 === 0 ? 1 : -1;
            const coefficient = this.binomialCoefficient(order, k);
            result.x += sign * coefficient * (A[order - k].x * Math.pow(w, k));
            result.y += sign * coefficient * (A[order - k].y * Math.pow(w, k));
            result.z += sign * coefficient * (A[order - k].z * Math.pow(w, k));
        }
        // Divide by w^(order+1)
        const wPower = Math.pow(w, order + 1);
        if (Math.abs(wPower) < Number.EPSILON) {
            throw new Error('Invalid derivative calculation: zero denominator');
        }
        const invWPower = 1 / wPower;
        const derivX = result.x * invWPower;
        const derivY = result.y * invWPower;
        const derivZ = result.z * invWPower;
        // Check for numerical validity
        if (!Number.isFinite(derivX) || !Number.isFinite(derivY) || !Number.isFinite(derivZ)) {
            throw new Error('Derivative calculation resulted in non-finite values');
        }
        // Return derivative vector
        return new Vector3D(derivX === 0 ? 0 : derivX, // Protect against -0
            derivY === 0 ? 0 : derivY, derivZ === 0 ? 0 : derivZ);
    }
    // ... other methods ...
    split(u) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Split parameter u must be a finite number');
        }
        if (u <= this.domain.min || u >= this.domain.max) {
            throw new Error('Split parameter must be within the curve domain (exclusive of endpoints)');
        }
        // Find the knot span containing u
        const knots = this.knotVector.knots;
        let insertionIndex = knots.findIndex(knot => knot > u);
        if (insertionIndex === -1) {
            throw new Error('Failed to find valid knot span for splitting');
        }
        // Calculate multiplicity of u in the knot vector
        let multiplicity = 0;
        for (let i = 0; i < knots.length; i++) {
            if (Math.abs(knots[i] - u) < Number.EPSILON) {
                multiplicity++;
            }
        }
        // Number of knots to insert (degree - multiplicity)
        const knotsToInsert = this.degree - multiplicity;
        // Insert knots at u to create full multiplicity
        let refinedKnots = [...knots];
        let refinedControlPoints = [...this.controlPoints];
        for (let insertion = 0; insertion < knotsToInsert; insertion++) {
            // Calculate alpha values for control point interpolation
            const alphas = [];
            for (let i = 0; i < this.degree - multiplicity; i++) {
                const span = insertionIndex - this.degree + i;
                if (span >= 0 && span < refinedKnots.length - 1) {
                    const denominator = refinedKnots[span + this.degree + 1] - refinedKnots[span];
                    const alpha = denominator > Number.EPSILON ? (u - refinedKnots[span]) / denominator : 0;
                    alphas.push(alpha);
                }
            }
            // Update control points
            const newControlPoints = [];
            for (let i = 0; i < insertionIndex - this.degree; i++) {
                newControlPoints.push(refinedControlPoints[i]);
            }
            for (let i = 0; i < this.degree - multiplicity + 1; i++) {
                const index = insertionIndex - this.degree + i;
                if (index >= 0 && index < refinedControlPoints.length - 1) {
                    const alpha = alphas[i];
                    const p1 = refinedControlPoints[index];
                    const p2 = refinedControlPoints[index + 1];
                    // Interpolate between control points
                    const x = (1 - alpha) * p1.x + alpha * p2.x;
                    const y = (1 - alpha) * p1.y + alpha * p2.y;
                    const z = (1 - alpha) * p1.z + alpha * p2.z;
                    const w = (1 - alpha) * p1.w + alpha * p2.w;
                    newControlPoints.push(new ControlPoint(x, y, z, w));
                }
            }
            for (let i = insertionIndex; i < refinedControlPoints.length; i++) {
                newControlPoints.push(refinedControlPoints[i]);
            }
            // Update knot vector
            const newKnots = [
                ...refinedKnots.slice(0, insertionIndex),
                u,
                ...refinedKnots.slice(insertionIndex)
            ];
            refinedKnots = newKnots;
            refinedControlPoints = newControlPoints;
            insertionIndex++;
        }
        // Create the two new curves by splitting at the full-multiplicity knot
        const splitIndex = refinedKnots.lastIndexOf(u);
        // Create first curve
        const firstKnots = refinedKnots.slice(0, splitIndex + this.degree + 1);
        const firstControlPoints = refinedControlPoints.slice(0, splitIndex);
        const firstCurve = new NurbsCurve(firstControlPoints, new KnotVector(firstKnots), this.degree);
        // Create second curve
        const secondKnots = refinedKnots.slice(splitIndex - this.degree);
        const secondControlPoints = refinedControlPoints.slice(splitIndex - this.degree);
        const secondCurve = new NurbsCurve(secondControlPoints, new KnotVector(secondKnots), this.degree);
        // Return both curves as an array
        return [
            firstCurve,
            secondCurve
        ];
    }
    // ... other methods ...
    join(otherCurve) {
        // Input validation
        if (!(otherCurve instanceof NurbsCurve)) {
            throw new Error('Join operation requires another NurbsCurve instance');
        }
        // Verify degrees match
        if (this.degree !== otherCurve.degree) {
            throw new Error('Cannot join curves of different degrees');
        }
        // Check if curves can be joined (end point of this curve matches start point of other curve)
        const thisEndPoint = this.evaluate(this.domain.max);
        const otherStartPoint = otherCurve.evaluate(otherCurve.domain.min);
        const joinPointTolerance = 1e-7;
        if (!thisEndPoint.equals(otherStartPoint, joinPointTolerance)) {
            throw new Error('Curves must share a common point for joining');
        }
        // Get control points from both curves
        const firstControlPoints = [...this.controlPoints];
        const secondControlPoints = [...otherCurve.controlPoints];
        // Get knot vectors from both curves
        const firstKnots = [...this.knotVector.knots];
        const secondKnots = [...otherCurve.knotVector.knots];
        // Scale and shift the second knot vector to align with the first
        const firstMax = this.domain.max;
        const secondMin = otherCurve.domain.min;
        const secondScale = otherCurve.domain.max - secondMin;
        const secondShift = firstMax - secondMin;
        const adjustedSecondKnots = secondKnots.map(knot => {
            const normalized = (knot - secondMin) / secondScale;
            return firstMax + normalized * secondScale;
        });
        // Merge knot vectors, removing duplicate knot at join point
        const mergedKnots = [
            ...firstKnots,
            ...adjustedSecondKnots.slice(this.degree + 1)
        ];
        // Merge control points, averaging coincident points at join
        const lastFirstIndex = firstControlPoints.length - 1;
        const firstSecondPoint = secondControlPoints[0];
        // Average the weights and positions of coincident control points
        const avgWeight = (firstControlPoints[lastFirstIndex].weight() + firstSecondPoint.weight()) / 2;
        const [x1, y1, z1] = firstControlPoints[lastFirstIndex].position();
        const [x2, y2, z2] = firstSecondPoint.position();
        const avgPoint = new ControlPoint((x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2, avgWeight);
        // Create merged control points array
        const mergedControlPoints = [
            ...firstControlPoints.slice(0, -1),
            avgPoint,
            ...secondControlPoints.slice(1)
        ];
        try {
            // Create new knot vector for joined curve
            const newKnotVector = new KnotVector(mergedKnots);
            // Validate the number of control points matches the knot vector
            const expectedControlPoints = newKnotVector.length - this.degree - 1;
            if (mergedControlPoints.length !== expectedControlPoints) {
                throw new Error('Invalid number of control points after joining');
            }
            // Create and return the joined curve
            return new NurbsCurve(mergedControlPoints, newKnotVector, this.degree);
        } catch (error) {
            throw new Error(`Failed to create joined curve: ${error.message}`);
        }
    }
    elevate() {
        // Create new degree
        const newDegree = this.degree + 1;
        // Get current control points in homogeneous form
        const currentPoints = this.controlPoints.map(cp => {
            const [x, y, z] = cp.position();
            const w = cp.weight();
            return [
                x * w,
                y * w,
                z * w,
                w
            ];
        });
        // Calculate new number of control points
        const n = this.controlPoints.length - 1;
        // current last index
        const newN = n + 1;
        // new last index
        // Initialize array for new control points
        const newPoints = new Array(newN + 1);
        // Calculate coefficient array
        const coefs = new Array(n + 1);
        for (let i = 0; i <= n; i++) {
            coefs[i] = this.binomialCoefficient(n, i) / this.binomialCoefficient(newN, i);
        }
        // Initialize first and last control points
        newPoints[0] = currentPoints[0];
        newPoints[newN] = currentPoints[n];
        // Calculate intermediate control points
        for (let i = 1; i < newN; i++) {
            // Initialize new point
            newPoints[i] = [
                0,
                0,
                0,
                0
            ];
            // Calculate alpha coefficients
            const alpha = i / newN;
            // Calculate new control point coordinates
            let denominator = 0;
            for (let j = Math.max(0, i - 1); j <= Math.min(n, i); j++) {
                const coef = coefs[j] * this.binomialCoefficient(i, j) * this.binomialCoefficient(newN - i, n - j) / this.binomialCoefficient(newN, n);
                // Add weighted contribution of current control point
                for (let k = 0; k < 4; k++) {
                    newPoints[i][k] += currentPoints[j][k] * coef;
                }
                denominator += coef;
            }
            // Normalize the point
            if (Math.abs(denominator) > Number.EPSILON) {
                for (let k = 0; k < 4; k++) {
                    newPoints[i][k] /= denominator;
                }
            }
        }
        // Convert homogeneous coordinates back to ControlPoint instances
        const elevatedControlPoints = newPoints.map(point => {
            const w = point[3];
            if (Math.abs(w) < Number.EPSILON) {
                throw new Error('Invalid weight in degree elevation calculation');
            }
            const invW = 1 / w;
            return new ControlPoint(point[0] * invW, point[1] * invW, point[2] * invW, w);
        });
        // Calculate new knot vector
        const newKnots = [];
        const oldKnots = this.knotVector.knots;
        // Copy first newDegree + 1 knots
        for (let i = 0; i <= newDegree; i++) {
            newKnots.push(oldKnots[0]);
        }
        // Copy middle knots
        for (let i = 1; i < oldKnots.length - this.degree - 1; i++) {
            newKnots.push(oldKnots[i]);
        }
        // Copy last newDegree + 1 knots
        for (let i = 0; i <= newDegree; i++) {
            newKnots.push(oldKnots[oldKnots.length - 1]);
        }
        try {
            // Create new knot vector
            const newKnotVector = new KnotVector(newKnots);
            // Return new elevated curve
            return new NurbsCurve(elevatedControlPoints, newKnotVector, newDegree);
        } catch (error) {
            throw new Error(`Failed to create elevated curve: ${error.message}`);
        }
    }
    // ... other methods ...
    reduce() {
        // Check if degree reduction is possible
        if (this.degree <= 1) {
            throw new Error('Cannot reduce curve of degree 1 or less');
        }
        const newDegree = this.degree - 1;
        const n = this.controlPoints.length - 1;
        // Current last index
        // Convert control points to homogeneous coordinates
        const currentPoints = this.controlPoints.map(cp => {
            const [x, y, z] = cp.position();
            const w = cp.weight();
            return [
                x * w,
                y * w,
                z * w,
                w
            ];
        });
        // Initialize matrices for the least squares problem
        const numNewPoints = n;
        const Q = new Array(numNewPoints + 1).fill(0).map(() => new Array(4).fill(0));
        const N = new Array(numNewPoints + 1).fill(0).map(() => new Array(numNewPoints + 1).fill(0));
        const R = new Array(numNewPoints + 1).fill(0).map(() => new Array(4).fill(0));
        // Sample parameters for least squares fitting
        const numSamples = Math.max(50, 5 * n);
        const parameterStep = (this.domain.max - this.domain.min) / (numSamples - 1);
        // Build least squares matrices
        for (let i = 0; i < numSamples; i++) {
            const u = this.domain.min + i * parameterStep;
            // Calculate basis functions for current and reduced degree
            const currentBasis = new Array(n + 1).fill(0);
            const reducedBasis = new Array(numNewPoints + 1).fill(0);
            // Calculate basis function values
            for (let j = 0; j <= n; j++) {
                currentBasis[j] = this.basisFunction.evaluate(j, this.degree, u);
            }
            for (let j = 0; j <= numNewPoints; j++) {
                reducedBasis[j] = this.basisFunction.evaluate(j, newDegree, u);
            }
            // Update Q matrix (control points times basis functions)
            for (let j = 0; j <= numNewPoints; j++) {
                for (let k = 0; k <= n; k++) {
                    const weight = currentBasis[k];
                    for (let d = 0; d < 4; d++) {
                        Q[j][d] += reducedBasis[j] * weight * currentPoints[k][d];
                    }
                }
            }
            // Update N matrix (basis functions products)
            for (let j = 0; j <= numNewPoints; j++) {
                for (let k = 0; k <= numNewPoints; k++) {
                    N[j][k] += reducedBasis[j] * reducedBasis[k];
                }
            }
        }
        // Solve the least squares system using LU decomposition
        try {
            // Solve NR = Q for each coordinate
            for (let d = 0; d < 4; d++) {
                const q = Q.map(row => row[d]);
                const r = this.solveLU(N, q);
                for (let i = 0; i <= numNewPoints; i++) {
                    R[i][d] = r[i];
                }
            }
        } catch (error) {
            throw new Error(`Failed to solve degree reduction system: ${error.message}`);
        }
        // Convert back from homogeneous coordinates to control points
        const reducedControlPoints = R.map(point => {
            const w = point[3];
            if (Math.abs(w) < Number.EPSILON) {
                throw new Error('Invalid weight in degree reduction calculation');
            }
            const invW = 1 / w;
            return new ControlPoint(point[0] * invW, point[1] * invW, point[2] * invW, w);
        });
        // Create new knot vector for reduced degree curve
        const oldKnots = this.knotVector.knots;
        const newKnots = [];
        // Copy first (newDegree + 1) knots
        for (let i = 0; i <= newDegree; i++) {
            newKnots.push(oldKnots[0]);
        }
        // Copy middle knots, removing every (degree - newDegree) knot
        for (let i = 1; i < oldKnots.length - this.degree - 1; i++) {
            newKnots.push(oldKnots[i]);
        }
        // Copy last (newDegree + 1) knots
        for (let i = 0; i <= newDegree; i++) {
            newKnots.push(oldKnots[oldKnots.length - 1]);
        }
        try {
            // Create new knot vector
            const newKnotVector = new KnotVector(newKnots);
            // Return new reduced curve
            return new NurbsCurve(reducedControlPoints, newKnotVector, newDegree);
        } catch (error) {
            throw new Error(`Failed to create reduced curve: ${error.message}`);
        }
    }
    // ... other methods ...
    reverse() {
        // Reverse control points
        const reversedControlPoints = [...this.controlPoints].reverse();
        // Get original knot vector values
        const oldKnots = this.knotVector.knots;
        // Calculate domain
        const a = oldKnots[0];
        const b = oldKnots[oldKnots.length - 1];
        // Create reversed knot vector
        const reversedKnots = new Float64Array(oldKnots.length);
        for (let i = 0; i < oldKnots.length; i++) {
            // Reverse and remap knots: t -> a + b - t
            reversedKnots[i] = a + b - oldKnots[oldKnots.length - 1 - i];
        }
        try {
            // Create new knot vector
            const newKnotVector = new KnotVector(reversedKnots);
            // Create and return the reversed curve
            return new NurbsCurve(reversedControlPoints, newKnotVector, this.degree);
        } catch (error) {
            throw new Error(`Failed to create reversed curve: ${error.message}`);
        }
    }
    closestPoint(point, tolerance = 1e-7, maxIterations = 100) {
        // Input validation
        if (!(point instanceof Vector3D)) {
            throw new Error('Input point must be a Vector3D instance');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        if (!Number.isInteger(maxIterations) || maxIterations <= 0) {
            throw new Error('Maximum iterations must be a positive integer');
        }
        // First do a coarse sampling to get initial parameter guess
        const numSamples = Math.max(20, this.controlPoints.length * 4);
        const parameterStep = (this.domain.max - this.domain.min) / (numSamples - 1);
        let minDistance = Infinity;
        let bestParameter = this.domain.min;
        // Sample points along the curve to find approximate closest point
        for (let i = 0; i < numSamples; i++) {
            const u = this.domain.min + i * parameterStep;
            const curvePoint = this.evaluate(u);
            const distance = point.subtract(curvePoint).lengthSquared();
            if (distance < minDistance) {
                minDistance = distance;
                bestParameter = u;
            }
        }
        // Refine the parameter using Newton's method
        let currentParameter = bestParameter;
        let converged = false;
        for (let iteration = 0; iteration < maxIterations; iteration++) {
            // Get point on curve and its derivatives
            const curvePoint = this.evaluate(currentParameter);
            const firstDerivative = this.derivative(currentParameter, 1);
            const secondDerivative = this.derivative(currentParameter, 2);
            // Calculate vector from curve point to target point
            const delta = point.subtract(curvePoint);
            // Calculate dot products for Newton iteration
            const dot1 = delta.dot(firstDerivative);
            const dot2 = firstDerivative.dot(firstDerivative) + delta.dot(secondDerivative);
            // Check if we're at a singular point
            if (Math.abs(dot2) < Number.EPSILON) {
                break;
            }
            // Calculate parameter update
            const step = dot1 / dot2;
            const newParameter = currentParameter + step;
            // Clamp parameter to domain
            const clampedParameter = Math.max(this.domain.min, Math.min(this.domain.max, newParameter));
            // Check for convergence
            if (Math.abs(clampedParameter - currentParameter) < tolerance) {
                converged = true;
                currentParameter = clampedParameter;
                break;
            }
            // Update parameter for next iteration
            currentParameter = clampedParameter;
        }
        // If Newton's method didn't converge, fall back to binary search
        if (!converged) {
            let left = this.domain.min;
            let right = this.domain.max;
            for (let i = 0; i < maxIterations && right - left > tolerance; i++) {
                const third = (right - left) / 3;
                const leftThird = left + third;
                const rightThird = right - third;
                const leftDistance = point.subtract(this.evaluate(leftThird)).lengthSquared();
                const rightDistance = point.subtract(this.evaluate(rightThird)).lengthSquared();
                if (leftDistance < rightDistance) {
                    right = rightThird;
                } else {
                    left = leftThird;
                }
            }
            currentParameter = (left + right) / 2;
        }
        // Return object containing parameter and point on curve
        return {
            parameter: currentParameter,
            point: this.evaluate(currentParameter),
            distance: point.subtract(this.evaluate(currentParameter)).length()
        };
    }
    // Helper method for calculating binomial coefficients
    binomialCoefficient(n, k) {
        if (k < 0 || k > n)
            return 0;
        if (k === 0 || k === n)
            return 1;
        let coefficient = 1;
        for (let i = 0; i < k; i++) {
            coefficient *= (n - i) / (i + 1);
        }
        return Math.round(coefficient);
    }
    // Helper method for solving linear system using LU decomposition
    solveLU(A, b) {
        const n = A.length;
        const L = new Array(n).fill(0).map(() => new Array(n).fill(0));
        const U = new Array(n).fill(0).map(() => new Array(n).fill(0));
        // LU decomposition
        for (let i = 0; i < n; i++) {
            // Upper triangular matrix
            for (let k = i; k < n; k++) {
                let sum = 0;
                for (let j = 0; j < i; j++) {
                    sum += L[i][j] * U[j][k];
                }
                U[i][k] = A[i][k] - sum;
            }
            // Lower triangular matrix
            for (let k = i; k < n; k++) {
                if (i === k) {
                    L[i][i] = 1;
                } else {
                    let sum = 0;
                    for (let j = 0; j < i; j++) {
                        sum += L[k][j] * U[j][i];
                    }
                    if (Math.abs(U[i][i]) < Number.EPSILON) {
                        throw new Error('Zero pivot encountered in LU decomposition');
                    }
                    L[k][i] = (A[k][i] - sum) / U[i][i];
                }
            }
        }
        // Forward substitution
        const y = new Array(n).fill(0);
        for (let i = 0; i < n; i++) {
            let sum = 0;
            for (let j = 0; j < i; j++) {
                sum += L[i][j] * y[j];
            }
            y[i] = b[i] - sum;
        }
        // Back substitution
        const x = new Array(n).fill(0);
        for (let i = n - 1; i >= 0; i--) {
            let sum = 0;
            for (let j = i + 1; j < n; j++) {
                sum += U[i][j] * x[j];
            }
            if (Math.abs(U[i][i]) < Number.EPSILON) {
                throw new Error('Zero pivot encountered in back substitution');
            }
            x[i] = (y[i] - sum) / U[i][i];
        }
        return x;
    }
}
class NurbsSurface {
    constructor() {
        this.controlPoints = [];
        this.knotVectorU = new KnotVector();
        this.knotVectorV = new KnotVector();
        this.degreeU = 0;
        this.degreeV = 0;
    }
    // ... other properties and methods ...
    evaluate(u, v) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (typeof v !== 'number' || !Number.isFinite(v)) {
            throw new Error('Parameter v must be a finite number');
        }
        // Check if parameters are within domain
        const uDomain = {
            min: this.knotVectorU.knots[this.degreeU],
            max: this.knotVectorU.knots[this.knotVectorU.length - this.degreeU - 1]
        };
        const vDomain = {
            min: this.knotVectorV.knots[this.degreeV],
            max: this.knotVectorV.knots[this.knotVectorV.length - this.degreeV - 1]
        };
        if (u < uDomain.min || u > uDomain.max) {
            throw new Error('Parameter u is outside the surface domain');
        }
        if (v < vDomain.min || v > vDomain.max) {
            throw new Error('Parameter v is outside the surface domain');
        }
        // Initialize basis function calculators
        const basisU = new BasisFunction();
        const basisV = new BasisFunction();
        basisU.knotVector = this.knotVectorU;
        basisV.knotVector = this.knotVectorV;
        basisU.degree = this.degreeU;
        basisV.degree = this.degreeV;
        // Initialize numerator and denominator for rational calculation
        let numeratorX = 0;
        let numeratorY = 0;
        let numeratorZ = 0;
        let denominator = 0;
        // Get surface dimensions
        const numU = this.controlPoints.length;
        const numV = this.controlPoints[0].length;
        // Calculate the sum using basis functions and control points
        for (let i = 0; i < numU; i++) {
            const Nu = basisU.evaluate(i, this.degreeU, u);
            for (let j = 0; j < numV; j++) {
                const Nv = basisV.evaluate(j, this.degreeV, v);
                // Get control point weight and coordinates
                const controlPoint = this.controlPoints[i][j];
                const weight = controlPoint.weight();
                const [x, y, z] = controlPoint.position();
                // Calculate weighted contribution
                const basisProduct = Nu * Nv * weight;
                // Add to numerator and denominator
                numeratorX += basisProduct * x;
                numeratorY += basisProduct * y;
                numeratorZ += basisProduct * z;
                denominator += basisProduct;
            }
        }
        // Check for zero denominator
        if (Math.abs(denominator) < Number.EPSILON) {
            throw new Error('Invalid surface evaluation: zero denominator');
        }
        // Calculate final point coordinates
        const invDenom = 1 / denominator;
        const x = numeratorX * invDenom;
        const y = numeratorY * invDenom;
        const z = numeratorZ * invDenom;
        // Check for numerical validity
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            throw new Error('Surface evaluation resulted in non-finite values');
        }
        // Return point on surface as Vector3D instance
        return new Vector3D(x === 0 ? 0 : x, // Protect against -0
            y === 0 ? 0 : y, z === 0 ? 0 : z);
    }
    derivative(u, v, orderU = 1, orderV = 1) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (typeof v !== 'number' || !Number.isFinite(v)) {
            throw new Error('Parameter v must be a finite number');
        }
        if (!Number.isInteger(orderU) || orderU < 0) {
            throw new Error('U derivative order must be a non-negative integer');
        }
        if (!Number.isInteger(orderV) || orderV < 0) {
            throw new Error('V derivative order must be a non-negative integer');
        }
        // If either order is greater than respective degree, return zero vector
        if (orderU > this.degreeU || orderV > this.degreeV) {
            return new Vector3D(0, 0, 0);
        }
        // Check parameter domains
        const uDomain = {
            min: this.knotVectorU.knots[this.degreeU],
            max: this.knotVectorU.knots[this.knotVectorU.length - this.degreeU - 1]
        };
        const vDomain = {
            min: this.knotVectorV.knots[this.degreeV],
            max: this.knotVectorV.knots[this.knotVectorV.length - this.degreeV - 1]
        };
        if (u < uDomain.min || u > uDomain.max) {
            throw new Error('Parameter u is outside the surface domain');
        }
        if (v < vDomain.min || v > vDomain.max) {
            throw new Error('Parameter v is outside the surface domain');
        }
        // Initialize basis function calculators
        const basisU = new BasisFunction();
        const basisV = new BasisFunction();
        basisU.knotVector = this.knotVectorU;
        basisV.knotVector = this.knotVectorV;
        basisU.degree = this.degreeU;
        basisV.degree = this.degreeV;
        // Get surface dimensions
        const numU = this.controlPoints.length;
        const numV = this.controlPoints[0].length;
        // Initialize arrays for numerator and denominator derivatives
        const derivs = {
            x: 0,
            y: 0,
            z: 0,
            w: 0
        };
        // Calculate basis function derivatives
        for (let i = 0; i < numU; i++) {
            const Nu = basisU.evaluateDerivative(i, this.degreeU, u, orderU);
            for (let j = 0; j < numV; j++) {
                const Nv = basisV.evaluateDerivative(j, this.degreeV, v, orderV);
                // Get control point data
                const controlPoint = this.controlPoints[i][j];
                const weight = controlPoint.weight();
                const [x, y, z] = controlPoint.position();
                // Calculate combined basis function derivative
                const basisDerivative = Nu * Nv * weight;
                // Add weighted contributions
                derivs.x += basisDerivative * x;
                derivs.y += basisDerivative * y;
                derivs.z += basisDerivative * z;
                derivs.w += basisDerivative;
            }
        }
        // Calculate point values for quotient rule
        const point = this.evaluate(u, v);
        const [px, py, pz] = [
            point.x,
            point.y,
            point.z
        ];
        // Apply quotient rule for rational surfaces
        // d(N/D) = (D*dN - N*dD) / D^2
        if (Math.abs(derivs.w) < Number.EPSILON) {
            throw new Error('Invalid derivative calculation: zero denominator');
        }
        const denomSqr = derivs.w * derivs.w;
        const x = (derivs.x * derivs.w - px * derivs.w * derivs.w) / denomSqr;
        const y = (derivs.y * derivs.w - py * derivs.w * derivs.w) / denomSqr;
        const z = (derivs.z * derivs.w - pz * derivs.w * derivs.w) / denomSqr;
        // Check for numerical validity
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            throw new Error('Derivative calculation resulted in non-finite values');
        }
        // Return derivative vector
        return new Vector3D(x === 0 ? 0 : x, // Protect against -0
            y === 0 ? 0 : y, z === 0 ? 0 : z);
    }
    // ... other methods ...
    normal(u, v) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (typeof v !== 'number' || !Number.isFinite(v)) {
            throw new Error('Parameter v must be a finite number');
        }
        try {
            // Calculate first partial derivatives
            const derivU = this.derivative(u, v, 1, 0);
            const derivV = this.derivative(u, v, 0, 1);
            // Check if either derivative is zero vector
            if (derivU.lengthSquared() < Number.EPSILON || derivV.lengthSquared() < Number.EPSILON) {
                throw new Error('Surface normal undefined: singular point');
            }
            // Calculate normal using cross product
            const normal = derivU.cross(derivV);
            // Check if normal is zero vector (indicates parallel derivatives)
            if (normal.lengthSquared() < Number.EPSILON) {
                throw new Error('Surface normal undefined: parallel derivatives');
            }
            // Return normalized vector
            return normal.normalize();
        } catch (error) {
            throw new Error(`Failed to calculate surface normal: ${error.message}`);
        }
    }
    // ... other methods ...
    split(parameter, direction = 'u') {
        // Input validation
        if (typeof parameter !== 'number' || !Number.isFinite(parameter)) {
            throw new Error('Split parameter must be a finite number');
        }
        if (direction !== 'u' && direction !== 'v') {
            throw new Error('Split direction must be either "u" or "v"');
        }
        // Get relevant knot vector and degree based on direction
        const knotVector = direction === 'u' ? this.knotVectorU : this.knotVectorV;
        const degree = direction === 'u' ? this.degreeU : this.degreeV;
        // Check if parameter is within domain
        const domain = {
            min: knotVector.knots[degree],
            max: knotVector.knots[knotVector.length - degree - 1]
        };
        if (parameter <= domain.min || parameter >= domain.max) {
            throw new Error('Split parameter must be within the surface domain (exclusive of endpoints)');
        }
        // Calculate multiplicity of parameter in knot vector
        let multiplicity = 0;
        const knots = knotVector.knots;
        for (const knot of knots) {
            if (Math.abs(knot - parameter) < Number.EPSILON) {
                multiplicity++;
            }
        }
        // Number of knots to insert to achieve full multiplicity
        const knotsToInsert = degree - multiplicity;
        if (knotsToInsert <= 0) {
            throw new Error('Surface is already split at this parameter');
        }
        // Perform knot insertion
        let refinedControlPoints = [...this.controlPoints];
        let refinedKnots = [...knots];
        for (let insertion = 0; insertion < knotsToInsert; insertion++) {
            // Find knot span
            let span = 0;
            while (span < refinedKnots.length - 1 && refinedKnots[span] <= parameter) {
                span++;
            }
            span--;
            // Create new arrays for the refined surface
            const newKnots = [
                ...refinedKnots.slice(0, span + 1),
                parameter,
                ...refinedKnots.slice(span + 1)
            ];
            const newControlPoints = direction === 'u' ? new Array(refinedControlPoints.length + 1).fill(null).map(() => []) : new Array(refinedControlPoints.length).fill(null).map(() => Array(refinedControlPoints[0].length + 1));
            // Calculate alpha values for control point interpolation
            const alphas = [];
            for (let i = 0; i < degree - multiplicity; i++) {
                const denominator = refinedKnots[span + degree + 1 - i] - refinedKnots[span + 1 - i];
                alphas[i] = denominator > Number.EPSILON ? (parameter - refinedKnots[span + 1 - i]) / denominator : 0;
            }
            // Update control points
            if (direction === 'u') {
                for (let j = 0; j < refinedControlPoints[0].length; j++) {
                    // Copy unaffected control points
                    for (let i = 0; i <= span - degree; i++) {
                        newControlPoints[i][j] = refinedControlPoints[i][j];
                    }
                    for (let i = span + 1; i < refinedControlPoints.length; i++) {
                        newControlPoints[i + 1][j] = refinedControlPoints[i][j];
                    }
                    // Calculate new control points
                    for (let i = 0; i <= degree - multiplicity; i++) {
                        const index = span - degree + i;
                        const alpha = alphas[i];
                        const p1 = refinedControlPoints[index][j];
                        const p2 = refinedControlPoints[index + 1][j];
                        // Interpolate between control points
                        const x = (1 - alpha) * p1.x + alpha * p2.x;
                        const y = (1 - alpha) * p1.y + alpha * p2.y;
                        const z = (1 - alpha) * p1.z + alpha * p2.z;
                        const w = (1 - alpha) * p1.w + alpha * p2.w;
                        newControlPoints[index + 1][j] = new ControlPoint(x, y, z, w);
                    }
                }
            } else {
                // v direction
                for (let i = 0; i < refinedControlPoints.length; i++) {
                    // Copy unaffected control points
                    for (let j = 0; j <= span - degree; j++) {
                        newControlPoints[i][j] = refinedControlPoints[i][j];
                    }
                    for (let j = span + 1; j < refinedControlPoints[0].length; j++) {
                        newControlPoints[i][j + 1] = refinedControlPoints[i][j];
                    }
                    // Calculate new control points
                    for (let j = 0; j <= degree - multiplicity; j++) {
                        const index = span - degree + j;
                        const alpha = alphas[j];
                        const p1 = refinedControlPoints[i][index];
                        const p2 = refinedControlPoints[i][index + 1];
                        // Interpolate between control points
                        const x = (1 - alpha) * p1.x + alpha * p2.x;
                        const y = (1 - alpha) * p1.y + alpha * p2.y;
                        const z = (1 - alpha) * p1.z + alpha * p2.z;
                        const w = (1 - alpha) * p1.w + alpha * p2.w;
                        newControlPoints[i][index + 1] = new ControlPoint(x, y, z, w);
                    }
                }
            }
            refinedControlPoints = newControlPoints;
            refinedKnots = newKnots;
        }
        // Create the two new surfaces by splitting at the parameter
        const splitIndex = refinedKnots.lastIndexOf(parameter);
        try {
            // Create first surface
            const firstSurface = new NurbsSurface();
            firstSurface.degreeU = this.degreeU;
            firstSurface.degreeV = this.degreeV;
            if (direction === 'u') {
                firstSurface.knotVectorU = new KnotVector(refinedKnots.slice(0, splitIndex + degree + 1));
                firstSurface.knotVectorV = this.knotVectorV;
                firstSurface.controlPoints = refinedControlPoints.slice(0, splitIndex);
            } else {
                firstSurface.knotVectorU = this.knotVectorU;
                firstSurface.knotVectorV = new KnotVector(refinedKnots.slice(0, splitIndex + degree + 1));
                firstSurface.controlPoints = refinedControlPoints.map(row => row.slice(0, splitIndex));
            }
            // Create second surface
            const secondSurface = new NurbsSurface();
            secondSurface.degreeU = this.degreeU;
            secondSurface.degreeV = this.degreeV;
            if (direction === 'u') {
                secondSurface.knotVectorU = new KnotVector(refinedKnots.slice(splitIndex - degree));
                secondSurface.knotVectorV = this.knotVectorV;
                secondSurface.controlPoints = refinedControlPoints.slice(splitIndex - degree);
            } else {
                secondSurface.knotVectorU = this.knotVectorU;
                secondSurface.knotVectorV = new KnotVector(refinedKnots.slice(splitIndex - degree));
                secondSurface.controlPoints = refinedControlPoints.map(row => row.slice(splitIndex - degree));
            }
            return [
                firstSurface,
                secondSurface
            ];
        } catch (error) {
            throw new Error(`Failed to create split surfaces: ${error.message}`);
        }
    }
    curvature(u, v) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (typeof v !== 'number' || !Number.isFinite(v)) {
            throw new Error('Parameter v must be a finite number');
        }
        try {
            // Calculate first and second derivatives
            const Su = this.derivative(u, v, 1, 0);
            const Sv = this.derivative(u, v, 0, 1);
            const Suu = this.derivative(u, v, 2, 0);
            const Svv = this.derivative(u, v, 0, 2);
            const Suv = this.derivative(u, v, 1, 1);
            // Calculate normal vector
            const normal = Su.cross(Sv).normalize();
            // Calculate coefficients of first fundamental form (E, F, G)
            const E = Su.dot(Su);
            const F = Su.dot(Sv);
            const G = Sv.dot(Sv);
            // Calculate coefficients of second fundamental form (L, M, N)
            const L = Suu.dot(normal);
            const M = Suv.dot(normal);
            const N = Svv.dot(normal);
            // Calculate determinants
            const firstFormDet = E * G - F * F;
            if (Math.abs(firstFormDet) < Number.EPSILON) {
                throw new Error('Degenerate surface point: zero first fundamental form determinant');
            }
            // Calculate shape operator matrix elements
            const shapeOp = {
                a11: (G * L - F * M) / firstFormDet,
                a12: (G * M - F * N) / firstFormDet,
                a21: (E * M - F * L) / firstFormDet,
                a22: (E * N - F * M) / firstFormDet
            };
            // Calculate eigenvalues (principal curvatures)
            const trace = shapeOp.a11 + shapeOp.a22;
            const det = shapeOp.a11 * shapeOp.a22 - shapeOp.a12 * shapeOp.a21;
            const discriminant = trace * trace - 4 * det;
            if (discriminant < -Number.EPSILON) {
                throw new Error('Complex principal curvatures encountered');
            }
            const sqrtDisc = Math.sqrt(Math.max(0, discriminant));
            const k1 = (trace + sqrtDisc) / 2;
            const k2 = (trace - sqrtDisc) / 2;
            // Calculate Gaussian and mean curvatures
            const K = det;
            // Gaussian curvature
            const H = trace / 2;
            // Mean curvature
            // Calculate principal directions if eigenvalues are distinct
            let dir1, dir2;
            if (Math.abs(k1 - k2) > Number.EPSILON) {
                // First principal direction
                const x1 = shapeOp.a12;
                const y1 = k1 - shapeOp.a11;
                dir1 = new Vector3D(x1, y1, 0).normalize();
                // Second principal direction
                const x2 = shapeOp.a12;
                const y2 = k2 - shapeOp.a11;
                dir2 = new Vector3D(x2, y2, 0).normalize();
            } else {
                // For umbilic points, any direction is principal
                dir1 = new Vector3D(1, 0, 0);
                dir2 = new Vector3D(0, 1, 0);
            }
            return {
                gaussianCurvature: K,
                meanCurvature: H,
                principalCurvatures: {
                    k1: k1,
                    k2: k2
                },
                principalDirections: {
                    dir1: dir1,
                    dir2: dir2
                },
                normal: normal
            };
        } catch (error) {
            throw new Error(`Failed to calculate surface curvature: ${error.message}`);
        }
    }
}
class TrimmedSurface {
    constructor() {
        this.baseSurface = null;
        this.trimCurves = [];
    }
    setTrimCurves(trimCurves, isOuter = true) {
        // Input validation
        if (!Array.isArray(trimCurves)) {
            throw new Error('Trim curves must be provided as an array');
        }
        // Validate each trim curve is a NURBS curve
        if (!trimCurves.every(curve => curve instanceof NurbsCurve)) {
            throw new Error('All trim curves must be instances of NurbsCurve');
        }
        // Validate base surface exists
        if (!this.baseSurface || !(this.baseSurface instanceof NurbsSurface)) {
            throw new Error('Base surface must be set before adding trim curves');
        }
        // Validate curve connectivity
        for (let i = 0; i < trimCurves.length; i++) {
            const currentCurve = trimCurves[i];
            const nextCurve = trimCurves[(i + 1) % trimCurves.length];
            // Get end points of current curve and start point of next curve
            const currentEnd = currentCurve.evaluate(currentCurve.domain.max);
            const nextStart = nextCurve.evaluate(nextCurve.domain.min);
            // Check if curves are connected (within tolerance)
            if (currentEnd.subtract(nextStart).length() > this.tolerance) {
                throw new Error(`Trim curves must form a closed loop: Gap detected between curves ${i} and ${(i + 1) % trimCurves.length}`);
            }
        }
        // Store trim curves with their orientation flag
        const trimLoop = {
            curves: trimCurves,
            isOuter: isOuter
        };
        // If this is the outer boundary, ensure it's the first in the list
        if (isOuter) {
            this.trimCurves.unshift(trimLoop);
        } else {
            this.trimCurves.push(trimLoop);
        }
        // Validate trim curve orientation
        this.validateTrimCurveOrientation(trimLoop);
        // Freeze the trim curve array to prevent direct modification
        Object.freeze(trimLoop.curves);
    }
    // ... other methods
    isPointInside(u, v) {
        // Input validation
        if (typeof u !== 'number' || !Number.isFinite(u)) {
            throw new Error('Parameter u must be a finite number');
        }
        if (typeof v !== 'number' || !Number.isFinite(v)) {
            throw new Error('Parameter v must be a finite number');
        }
        // Validate base surface exists
        if (!this.baseSurface || !(this.baseSurface instanceof NurbsSurface)) {
            throw new Error('Base surface must be set before testing point containment');
        }
        // Check if point is within surface parameter domain
        const uDomain = {
            min: this.baseSurface.knotVectorU.knots[this.baseSurface.degreeU],
            max: this.baseSurface.knotVectorU.knots[this.baseSurface.knotVectorU.length - this.baseSurface.degreeU - 1]
        };
        const vDomain = {
            min: this.baseSurface.knotVectorV.knots[this.baseSurface.degreeV],
            max: this.baseSurface.knotVectorV.knots[this.baseSurface.knotVectorV.length - this.baseSurface.degreeV - 1]
        };
        if (u < uDomain.min || u > uDomain.max || v < vDomain.min || v > vDomain.max) {
            return false;
        }
        // If no trim curves, point is inside if it's within domain
        if (this.trimCurves.length === 0) {
            return true;
        }
        // Test point in parameter space
        const testPoint = [
            u,
            v
        ];
        let isInside = false;
        // Process each trim loop
        for (const trimLoop of this.trimCurves) {
            let windingNumber = 0;
            const curves = trimLoop.curves;
            // Process each curve in the trim loop
            for (const curve of curves) {
                const numSamples = 50;
                // Number of samples for curve approximation
                const parameterStep = (curve.domain.max - curve.domain.min) / (numSamples - 1);
                // Process each segment of the curve
                for (let i = 0; i < numSamples - 1; i++) {
                    const t1 = curve.domain.min + i * parameterStep;
                    const t2 = curve.domain.min + (i + 1) * parameterStep;
                    const p1 = curve.evaluate(t1);
                    const p2 = curve.evaluate(t2);
                    // Convert to parameter space points
                    const p1Param = [
                        p1.x,
                        p1.y
                    ];
                    const p2Param = [
                        p2.x,
                        p2.y
                    ];
                    // Calculate winding number contribution
                    if (p1Param[1] <= testPoint[1]) {
                        if (p2Param[1] > testPoint[1] && this.isLeft(p1Param, p2Param, testPoint) > 0) {
                            windingNumber++;
                        }
                    } else {
                        if (p2Param[1] <= testPoint[1] && this.isLeft(p1Param, p2Param, testPoint) < 0) {
                            windingNumber--;
                        }
                    }
                }
            }
            // Update inside/outside status based on winding number and boundary type
            if (trimLoop.isOuter) {
                // For outer boundary, non-zero winding number means inside
                isInside = windingNumber !== 0;
            } else {
                // For inner boundaries (holes), non-zero winding number means outside the hole
                if (windingNumber !== 0) {
                    return false;
                }
            }
        }
        return isInside;
    }
    createHole(trimCurves) {
        // Input validation
        if (!Array.isArray(trimCurves)) {
            throw new Error('Trim curves must be provided as an array');
        }
        // Verify all curves are NURBS curves
        if (!trimCurves.every(curve => curve instanceof NurbsCurve)) {
            throw new Error('All trim curves must be instances of NurbsCurve');
        }
        // Validate base surface exists
        if (!this.baseSurface || !(this.baseSurface instanceof NurbsSurface)) {
            throw new Error('Base surface must be set before creating holes');
        }
        // Verify the trim curves form a closed loop
        for (let i = 0; i < trimCurves.length; i++) {
            const currentCurve = trimCurves[i];
            const nextCurve = trimCurves[(i + 1) % trimCurves.length];
            // Get end points
            const currentEnd = currentCurve.evaluate(currentCurve.domain.max);
            const nextStart = nextCurve.evaluate(nextCurve.domain.min);
            // Check connectivity within tolerance
            if (currentEnd.subtract(nextStart).length() > this.tolerance) {
                throw new Error(`Trim curves must form a closed loop: Gap detected between curves ${i} and ${(i + 1) % trimCurves.length}`);
            }
        }
        for (const existingLoop of this.trimCurves) {
            for (const existingCurve of existingLoop.curves) {
                for (const newCurve of trimCurves) {
                    // Create intersection calculator
                    const intersectionCalc = new IntersectionCalculator();
                    try {
                        const intersections = intersectionCalc.findCurveCurveIntersections(existingCurve, newCurve);
                        if (intersections && intersections.points.length > 0) {
                            throw new Error('New hole intersects with existing trim curves');
                        }
                    } catch (error) {
                        throw new Error(`Failed to check curve intersections: ${error.message}`);
                    }
                }
            }
        }
        // Verify the hole is within the outer boundary (if one exists)
        if (this.trimCurves.length > 0 && this.trimCurves[0].isOuter) {
            // Sample points along the new hole boundary
            const numSamplePoints = 20;
            for (const curve of trimCurves) {
                const parameterStep = (curve.domain.max - curve.domain.min) / numSamplePoints;
                for (let i = 0; i <= numSamplePoints; i++) {
                    const parameter = curve.domain.min + i * parameterStep;
                    const point = curve.evaluate(parameter);
                    // Check if point is inside the outer boundary
                    const isInsideOuter = this.isPointInside(point.x, point.y);
                    if (!isInsideOuter) {
                        throw new Error('Hole must be completely within the outer boundary');
                    }
                }
            }
        }
        // Create the new trim loop for the hole
        const holeLoop = {
            curves: Object.freeze([...trimCurves]),
            isOuter: false
        };
        // Validate and correct the orientation of the hole's trim curves
        this.validateTrimCurveOrientation(holeLoop);
        // Add the hole to the trim curves array
        this.trimCurves.push(holeLoop);
        // Return the index of the newly created hole
        return this.trimCurves.length - 1;
    }
    tessellate(surfaceTolerance = 0.001, trimTolerance = 0.0001) {
        // Input validation
        if (!this.baseSurface || !(this.baseSurface instanceof NurbsSurface)) {
            throw new Error('Base surface must be set before tessellation');
        }
        if (typeof surfaceTolerance !== 'number' || surfaceTolerance <= 0) {
            throw new Error('Surface tolerance must be a positive number');
        }
        if (typeof trimTolerance !== 'number' || trimTolerance <= 0) {
            throw new Error('Trim tolerance must be a positive number');
        }
        // Create initial surface tessellation grid
        const tessellator = new Tessellator();
        const uDomain = {
            min: this.baseSurface.knotVectorU.knots[this.baseSurface.degreeU],
            max: this.baseSurface.knotVectorU.knots[this.baseSurface.knotVectorU.length - this.baseSurface.degreeU - 1]
        };
        const vDomain = {
            min: this.baseSurface.knotVectorV.knots[this.baseSurface.degreeV],
            max: this.baseSurface.knotVectorV.knots[this.baseSurface.knotVectorV.length - this.baseSurface.degreeV - 1]
        };
        // Calculate initial grid size based on surface size and curvature
        const initialGridSize = this.calculateInitialGridSize(uDomain, vDomain);
        const uStep = (uDomain.max - uDomain.min) / initialGridSize.u;
        const vStep = (vDomain.max - vDomain.min) / initialGridSize.v;
        // Generate initial grid points
        const gridPoints = [];
        const gridParams = [];
        for (let i = 0; i <= initialGridSize.u; i++) {
            const row = [];
            const paramRow = [];
            const u = uDomain.min + i * uStep;
            for (let j = 0; j <= initialGridSize.v; j++) {
                const v = vDomain.min + j * vStep;
                if (this.isPointInside(u, v)) {
                    const point = this.baseSurface.evaluate(u, v);
                    row.push(point);
                    paramRow.push({
                        u,
                        v
                    });
                } else {
                    row.push(null);
                    paramRow.push(null);
                }
            }
            gridPoints.push(row);
            gridParams.push(paramRow);
        }
        // Refine grid based on surface curvature
        this.refineGrid(gridPoints, gridParams, surfaceTolerance);
        // Process trim curves
        const trimCurvePoints = [];
        for (const trimLoop of this.trimCurves) {
            const loopPoints = [];
            for (const curve of trimLoop.curves) {
                const tessResult = tessellator.tessellateNurbsCurve(curve, trimTolerance);
                loopPoints.push({
                    points: tessResult.points,
                    parameters: tessResult.parameters,
                    isOuter: trimLoop.isOuter
                });
            }
            trimCurvePoints.push(loopPoints);
        }
        // Generate triangles from grid points
        const triangles = this.generateTriangles(gridPoints, gridParams);
        // Return tessellation result
        return {
            points: this.flattenPoints(gridPoints),
            parameters: this.flattenPoints(gridParams),
            triangles: triangles,
            trimCurves: trimCurvePoints
        };
    }
    validateTrimCurveOrientation(trimLoop) {
        // Calculate approximate signed area of the trim loop
        let signedArea = 0;
        const curves = trimLoop.curves;
        for (let i = 0; i < curves.length; i++) {
            const curve = curves[i];
            const numSamples = 20;
            // Number of samples per curve
            const parameterStep = (curve.domain.max - curve.domain.min) / (numSamples - 1);
            for (let j = 0; j < numSamples - 1; j++) {
                const t1 = curve.domain.min + j * parameterStep;
                const t2 = curve.domain.min + (j + 1) * parameterStep;
                const p1 = curve.evaluate(t1);
                const p2 = curve.evaluate(t2);
                // Use Green's theorem to calculate signed area
                signedArea += (p2.x - p1.x) * (p2.y + p1.y) / 2;
            }
        }
        // For outer boundaries, signed area should be positive
        // For inner boundaries (holes), signed area should be negative
        const correctOrientation = trimLoop.isOuter ? signedArea > 0 : signedArea < 0;
        if (!correctOrientation) {
            // Reverse the curves if orientation is incorrect
            trimLoop.curves.reverse();
            for (let i = 0; i < trimLoop.curves.length; i++) {
                trimLoop.curves[i] = trimLoop.curves[i].reverse();
            }
        }
    }
    // Helper method to determine if test point is left/right/on the line from p1 to p2
    isLeft(p1, p2, testPoint) {
        return (p2[0] - p1[0]) * (testPoint[1] - p1[1]) - (testPoint[0] - p1[0]) * (p2[1] - p1[1]);
    }
    calculateInitialGridSize(uDomain, vDomain) {
        // Calculate initial grid size based on surface properties
        const uSpan = uDomain.max - uDomain.min;
        const vSpan = vDomain.max - vDomain.min;
        // Sample surface curvature at several points
        const numSamples = 9;
        let maxCurvature = 0;
        for (let i = 0; i < numSamples; i++) {
            const u = uDomain.min + uSpan * i / (numSamples - 1);
            for (let j = 0; j < numSamples; j++) {
                const v = vDomain.min + vSpan * j / (numSamples - 1);
                try {
                    const curvatureInfo = this.baseSurface.curvature(u, v);
                    const k1 = Math.abs(curvatureInfo.principalCurvatures.k1);
                    const k2 = Math.abs(curvatureInfo.principalCurvatures.k2);
                    maxCurvature = Math.max(maxCurvature, k1, k2);
                } catch (error) {
                    // Skip problematic curvature calculations
                    continue;
                }
            }
        }
        // Calculate grid size based on curvature and span
        const baseSize = 10;
        // Minimum grid size
        const curvatureFactor = Math.sqrt(maxCurvature) * 100;
        return {
            u: Math.max(baseSize, Math.ceil(uSpan * curvatureFactor)),
            v: Math.max(baseSize, Math.ceil(vSpan * curvatureFactor))
        };
    }
    refineGrid(gridPoints, gridParams, tolerance) {
        // Refine grid based on surface properties
        let needsRefinement;
        do {
            needsRefinement = false;
            // Check each grid cell for refinement
            for (let i = 0; i < gridPoints.length - 1; i++) {
                for (let j = 0; j < gridPoints[i].length - 1; j++) {
                    if (this.needsRefinement(gridPoints[i][j], gridPoints[i + 1][j], gridPoints[i][j + 1], gridPoints[i + 1][j + 1], gridParams[i][j], gridParams[i + 1][j], gridParams[i][j + 1], gridParams[i + 1][j + 1], tolerance)) {
                        this.subdivideCell(i, j, gridPoints, gridParams);
                        needsRefinement = true;
                    }
                }
            }
        } while (needsRefinement);
    }
    needsRefinement(p1, p2, p3, p4, param1, param2, param3, param4, tolerance) {
        if (!p1 || !p2 || !p3 || !p4)
            return false;
        // Calculate midpoint parameters
        const midU = (param1.u + param2.u) / 2;
        const midV = (param1.v + param3.v) / 2;
        // Evaluate surface at midpoint
        const actualMidpoint = this.baseSurface.evaluate(midU, midV);
        // Calculate interpolated midpoint
        const interpolatedMidpoint = p1.add(p2).add(p3).add(p4).multiply(0.25);
        // Check deviation
        return actualMidpoint.subtract(interpolatedMidpoint).length() > tolerance;
    }
    generateTriangles(gridPoints, gridParams) {
        const triangles = [];
        for (let i = 0; i < gridPoints.length - 1; i++) {
            for (let j = 0; j < gridPoints[i].length - 1; j++) {
                const p00 = gridPoints[i][j];
                const p10 = gridPoints[i + 1][j];
                const p01 = gridPoints[i][j + 1];
                const p11 = gridPoints[i + 1][j + 1];
                // Only create triangles if all points exist
                if (p00 && p10 && p01 && p11) {
                    // Add two triangles for the quad
                    triangles.push([
                        this.getPointIndex(p00, gridPoints),
                        this.getPointIndex(p10, gridPoints),
                        this.getPointIndex(p01, gridPoints)
                    ]);
                    triangles.push([
                        this.getPointIndex(p10, gridPoints),
                        this.getPointIndex(p11, gridPoints),
                        this.getPointIndex(p01, gridPoints)
                    ]);
                }
            }
        }
        return triangles;
    }
    flattenPoints(grid) {
        const flattened = [];
        for (const row of grid) {
            for (const point of row) {
                if (point !== null) {
                    flattened.push(point);
                }
            }
        }
        return flattened;
    }
    getPointIndex(point, gridPoints) {
        let index = 0;
        for (const row of gridPoints) {
            for (const gridPoint of row) {
                if (gridPoint !== null) {
                    if (gridPoint === point) {
                        return index;
                    }
                    index++;
                }
            }
        }
        return -1;
    }
}
class IntersectionCalculator {
    constructor() {
        this.tolerance = defaultTolerance;
    }
    findCurveCurveIntersections(curve1, curve2, tolerance = this.tolerance) {
        // Input validation
        if (!(curve1 instanceof NurbsCurve) || !(curve2 instanceof NurbsCurve)) {
            throw new Error('Both inputs must be NURBS curves');
        }
        // Initialize result object
        const result = new IntersectionResult();
        result.points = [];
        result.parameters = [];
        // Check for overlapping bounding boxes first
        const box1 = this.computeBoundingBox(curve1);
        const box2 = this.computeBoundingBox(curve2);
        if (!this.boundingBoxesIntersect(box1, box2)) {
            return result;
        }
        // Initialize subdivision intervals
        const intervals1 = [{
            t1: curve1.domain.min,
            t2: curve1.domain.max
        }];
        const intervals2 = [{
            t1: curve2.domain.min,
            t2: curve2.domain.max
        }];
        // Process subdivided intervals
        const candidateIntervals = [];
        while (intervals1.length > 0 && intervals2.length > 0) {
            const i1 = intervals1.pop();
            const i2 = intervals2.pop();
            // Compute bounding boxes for current intervals
            const box1 = this.computeIntervalBoundingBox(curve1, i1.t1, i1.t2);
            const box2 = this.computeIntervalBoundingBox(curve2, i2.t1, i2.t2);
            // Check if boxes intersect
            if (this.boundingBoxesIntersect(box1, box2)) {
                // Check interval sizes
                const size1 = i1.t2 - i1.t1;
                const size2 = i2.t2 - i2.t1;
                if (size1 < tolerance && size2 < tolerance) {
                    // Intervals are small enough - add as candidate intersection
                    candidateIntervals.push({
                        t1: (i1.t1 + i1.t2) / 2,
                        t2: (i2.t1 + i2.t2) / 2
                    });
                } else {
                    // Subdivide larger interval
                    if (size1 > size2) {
                        const tmid = (i1.t1 + i1.t2) / 2;
                        intervals1.push({
                            t1: i1.t1,
                            t2: tmid
                        }, {
                            t1: tmid,
                            t2: i1.t2
                        });
                        intervals2.push(i2, i2);
                    } else {
                        const tmid = (i2.t1 + i2.t2) / 2;
                        intervals2.push({
                            t1: i2.t1,
                            t2: tmid
                        }, {
                            t1: tmid,
                            t2: i2.t2
                        });
                        intervals1.push(i1, i1);
                    }
                }
            }
        }
        // Refine candidate intersections using Newton iteration
        const processedPoints = new Set();
        for (const candidate of candidateIntervals) {
            try {
                const refined = this.refineIntersection(curve1, curve2, candidate.t1, candidate.t2, tolerance);
                if (refined) {
                    const point = curve1.evaluate(refined.t1);
                    const pointKey = `${point.x.toFixed(6)},${point.y.toFixed(6)},${point.z.toFixed(6)}`;
                    // Check if this point is already found
                    if (!processedPoints.has(pointKey)) {
                        processedPoints.add(pointKey);
                        result.points.push(point);
                        result.parameters.push({
                            curve1: refined.t1,
                            curve2: refined.t2
                        });
                    }
                }
            } catch (error) {
                // Skip failed refinements
                continue;
            }
        }
        // Sort intersections by curve1 parameter
        const sorted = result.parameters.map((param, index) => ({
            param,
            index
        })).sort((a, b) => a.param.curve1 - b.param.curve1);
        result.points = sorted.map(item => result.points[item.index]);
        result.parameters = sorted.map(item => result.parameters[item.index]);
        return result;
    }
    findSurfaceSurfaceIntersections(surface1, surface2, tolerance = this.tolerance) {
        // Input validation
        if (!(surface1 instanceof NurbsSurface) || !(surface2 instanceof NurbsSurface)) {
            throw new Error('Both inputs must be NURBS surfaces');
        }
        // Initialize result object
        const result = new IntersectionResult();
        result.points = [];
        result.curves = [];
        result.parameters = [];
        // Check for overlapping bounding boxes first
        const box1 = this.computeBoundingBox(surface1);
        const box2 = this.computeBoundingBox(surface2);
        if (!this.boundingBoxesIntersect(box1, box2)) {
            return result;
        }
        // Initial subdivision of parameter spaces
        const subdivisions1 = [{
            u1: surface1.knotVectorU.knots[surface1.degreeU],
            u2: surface1.knotVectorU.knots[surface1.knotVectorU.length - surface1.degreeU - 1],
            v1: surface1.knotVectorV.knots[surface1.degreeV],
            v2: surface1.knotVectorV.knots[surface1.knotVectorV.length - surface1.degreeV - 1]
        }];
        const subdivisions2 = [{
            u1: surface2.knotVectorU.knots[surface2.degreeU],
            u2: surface2.knotVectorU.knots[surface2.knotVectorU.length - surface2.degreeU - 1],
            v1: surface2.knotVectorV.knots[surface2.degreeV],
            v2: surface2.knotVectorV.knots[surface2.knotVectorV.length - surface2.degreeV - 1]
        }];
        // Arrays to store intersection points and curve segments
        const intersectionPoints = new Set();
        const curveSegments = [];
        // Process subdivided patches
        while (subdivisions1.length > 0 && subdivisions2.length > 0) {
            const patch1 = subdivisions1.pop();
            const patch2 = subdivisions2.pop();
            // Compute bounding boxes for current patches
            const patchBox1 = this.computePatchBoundingBox(surface1, patch1);
            const patchBox2 = this.computePatchBoundingBox(surface2, patch2);
            if (this.boundingBoxesIntersect(patchBox1, patchBox2)) {
                // Calculate patch sizes
                const size1 = Math.max(patch1.u2 - patch1.u1, patch1.v2 - patch1.v1);
                const size2 = Math.max(patch2.u2 - patch2.u1, patch2.v2 - patch2.v1);
                if (size1 < tolerance && size2 < tolerance) {
                    // Patches are small enough - compute intersection point
                    const p1 = surface1.evaluate((patch1.u1 + patch1.u2) / 2, (patch1.v1 + patch1.v2) / 2);
                    const p2 = surface2.evaluate((patch2.u1 + patch2.u2) / 2, (patch2.v1 + patch2.v2) / 2);
                    if (p1.subtract(p2).length() < tolerance) {
                        // Found an intersection point
                        const params = {
                            surface1: {
                                u: (patch1.u1 + patch1.u2) / 2,
                                v: (patch1.v1 + patch1.v2) / 2
                            },
                            surface2: {
                                u: (patch2.u1 + patch2.u2) / 2,
                                v: (patch2.v1 + patch2.v2) / 2
                            }
                        };
                        const pointKey = `${p1.x.toFixed(6)},${p1.y.toFixed(6)},${p1.z.toFixed(6)}`;
                        if (!intersectionPoints.has(pointKey)) {
                            intersectionPoints.add(pointKey);
                            result.points.push(p1);
                            result.parameters.push(params);
                            curveSegments.push({
                                point: p1,
                                params: params
                            });
                        }
                    }
                } else {
                    // Subdivide larger patch
                    if (size1 > size2) {
                        // Subdivide patch1
                        const uMid = (patch1.u1 + patch1.u2) / 2;
                        const vMid = (patch1.v1 + patch1.v2) / 2;
                        subdivisions1.push({
                            u1: patch1.u1,
                            u2: uMid,
                            v1: patch1.v1,
                            v2: vMid
                        }, {
                            u1: uMid,
                            u2: patch1.u2,
                            v1: patch1.v1,
                            v2: vMid
                        }, {
                            u1: patch1.u1,
                            u2: uMid,
                            v1: vMid,
                            v2: patch1.v2
                        }, {
                            u1: uMid,
                            u2: patch1.u2,
                            v1: vMid,
                            v2: patch1.v2
                        });
                        subdivisions2.push(patch2, patch2, patch2, patch2);
                    } else {
                        // Subdivide patch2
                        const uMid = (patch2.u1 + patch2.u2) / 2;
                        const vMid = (patch2.v1 + patch2.v2) / 2;
                        subdivisions2.push({
                            u1: patch2.u1,
                            u2: uMid,
                            v1: patch2.v1,
                            v2: vMid
                        }, {
                            u1: uMid,
                            u2: patch2.u2,
                            v1: patch2.v1,
                            v2: vMid
                        }, {
                            u1: patch2.u1,
                            u2: uMid,
                            v1: vMid,
                            v2: patch2.v2
                        }, {
                            u1: uMid,
                            u2: patch2.u2,
                            v1: vMid,
                            v2: patch2.v2
                        });
                        subdivisions1.push(patch1, patch1, patch1, patch1);
                    }
                }
            }
        }
        // Create intersection curves from segments
        if (curveSegments.length > 0) {
            const curves = this.reconstructIntersectionCurves(curveSegments, tolerance);
            result.curves = curves;
        }
        return result;
    }
    // ... other methods ...
    findSurfaceCurveIntersections(surface, curve, tolerance = this.tolerance) {
        // Input validation
        if (!(surface instanceof NurbsSurface)) {
            throw new Error('First parameter must be a NURBS surface');
        }
        if (!(curve instanceof NurbsCurve)) {
            throw new Error('Second parameter must be a NURBS curve');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        // Initialize result object
        const result = new IntersectionResult();
        result.points = [];
        result.parameters = [];
        // Check bounding boxes first
        const surfaceBox = this.computeBoundingBox(surface);
        const curveBox = this.computeBoundingBox(curve);
        if (!this.boundingBoxesIntersect(surfaceBox, curveBox)) {
            return result;
        }
        // Initialize subdivision intervals for the curve
        const curveIntervals = [{
            t1: curve.domain.min,
            t2: curve.domain.max
        }];
        // Initialize subdivision patches for the surface
        const surfacePatches = [{
            u1: surface.knotVectorU.knots[surface.degreeU],
            u2: surface.knotVectorU.knots[surface.knotVectorU.length - surface.degreeU - 1],
            v1: surface.knotVectorV.knots[surface.degreeV],
            v2: surface.knotVectorV.knots[surface.knotVectorV.length - surface.degreeV - 1]
        }];
        // Store candidate intersection points
        const candidates = [];
        // Process subdivided intervals and patches
        while (curveIntervals.length > 0 && surfacePatches.length > 0) {
            const interval = curveIntervals.pop();
            const patch = surfacePatches.pop();
            // Compute bounding boxes for current interval and patch
            const curveSegmentBox = this.computeIntervalBoundingBox(curve, interval.t1, interval.t2);
            const surfacePatchBox = this.computePatchBoundingBox(surface, patch);
            if (this.boundingBoxesIntersect(curveSegmentBox, surfacePatchBox)) {
                // Calculate sizes of curve segment and surface patch
                const curveSize = interval.t2 - interval.t1;
                const surfaceSize = Math.max(patch.u2 - patch.u1, patch.v2 - patch.v1);
                if (curveSize < tolerance && surfaceSize < tolerance) {
                    // Compute candidate intersection point
                    const t = (interval.t1 + interval.t2) / 2;
                    const u = (patch.u1 + patch.u2) / 2;
                    const v = (patch.v1 + patch.v2) / 2;
                    const curvePoint = curve.evaluate(t);
                    const surfacePoint = surface.evaluate(u, v);
                    if (curvePoint.subtract(surfacePoint).length() < tolerance) {
                        candidates.push({
                            point: curvePoint,
                            curveParam: t,
                            surfaceParams: {
                                u,
                                v
                            }
                        });
                    }
                } else {
                    // Subdivide the larger of curve or surface
                    if (curveSize > surfaceSize) {
                        // Subdivide curve interval
                        const tmid = (interval.t1 + interval.t2) / 2;
                        curveIntervals.push({
                            t1: interval.t1,
                            t2: tmid
                        }, {
                            t1: tmid,
                            t2: interval.t2
                        });
                        surfacePatches.push(patch, patch);
                    } else {
                        // Subdivide surface patch
                        const umid = (patch.u1 + patch.u2) / 2;
                        const vmid = (patch.v1 + patch.v2) / 2;
                        surfacePatches.push({
                            u1: patch.u1,
                            u2: umid,
                            v1: patch.v1,
                            v2: vmid
                        }, {
                            u1: umid,
                            u2: patch.u2,
                            v1: patch.v1,
                            v2: vmid
                        }, {
                            u1: patch.u1,
                            u2: umid,
                            v1: vmid,
                            v2: patch.v2
                        }, {
                            u1: umid,
                            u2: patch.u2,
                            v1: vmid,
                            v2: patch.v2
                        });
                        curveIntervals.push(interval, interval, interval, interval);
                    }
                }
            }
        }
        // Refine candidate intersections using Newton iteration
        const processedPoints = new Set();
        for (const candidate of candidates) {
            try {
                const refined = this.refineSurfaceCurveIntersection(surface, curve, candidate.surfaceParams.u, candidate.surfaceParams.v, candidate.curveParam, tolerance);
                if (refined) {
                    const point = curve.evaluate(refined.t);
                    const pointKey = `${point.x.toFixed(6)},${point.y.toFixed(6)},${point.z.toFixed(6)}`;
                    if (!processedPoints.has(pointKey)) {
                        processedPoints.add(pointKey);
                        result.points.push(point);
                        result.parameters.push({
                            curve: refined.t,
                            surface: {
                                u: refined.u,
                                v: refined.v
                            }
                        });
                    }
                }
            } catch (error) {
                // Skip failed refinements
                continue;
            }
        }
        // Sort intersections by curve parameter
        const sorted = result.parameters.map((param, index) => ({
            param,
            index
        })).sort((a, b) => a.param.curve - b.param.curve);
        result.points = sorted.map(item => result.points[item.index]);
        result.parameters = sorted.map(item => result.parameters[item.index]);
        return result;
    }
    refinementIntersection(curve1, curve2, t1, t2, tolerance = this.tolerance, maxIterations = 20) {
        // Input validation
        if (!(curve1 instanceof NurbsCurve) || !(curve2 instanceof NurbsCurve)) {
            throw new Error('Both inputs must be NURBS curves');
        }
        if (typeof t1 !== 'number' || typeof t2 !== 'number') {
            throw new Error('Parameters t1 and t2 must be numbers');
        }
        if (!Number.isFinite(t1) || !Number.isFinite(t2)) {
            throw new Error('Parameters t1 and t2 must be finite numbers');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        if (!Number.isInteger(maxIterations) || maxIterations <= 0) {
            throw new Error('Maximum iterations must be a positive integer');
        }
        let currentT1 = t1;
        let currentT2 = t2;
        try {
            // Newton-Raphson iteration
            for (let iteration = 0; iteration < maxIterations; iteration++) {
                // Get points and derivatives at current parameters
                const p1 = curve1.evaluate(currentT1);
                const p2 = curve2.evaluate(currentT2);
                const d1 = curve1.derivative(currentT1, 1);
                const d2 = curve2.derivative(currentT2, 1);
                // Calculate distance vector between points
                const distance = p2.subtract(p1);
                // Check if points are close enough
                if (distance.length() < tolerance) {
                    return {
                        t1: currentT1,
                        t2: currentT2,
                        point: p1,
                        distance: distance.length()
                    };
                }
                // Set up system of equations for Newton iteration
                // [d1.x  -d2.x] [dt1] = [distance.x]
                // [d1.y  -d2.y] [dt2] = [distance.y]
                // [d1.z  -d2.z]        [distance.z]
                // Use least squares to solve overdetermined system
                const A = [
                    [
                        d1.x,
                        -d2.x
                    ],
                    [
                        d1.y,
                        -d2.y
                    ],
                    [
                        d1.z,
                        -d2.z
                    ]
                ];
                const b = [
                    distance.x,
                    distance.y,
                    distance.z
                ];
                // Calculate normal equations
                const AtA = [
                    [
                        A[0][0] * A[0][0] + A[1][0] * A[1][0] + A[2][0] * A[2][0],
                        A[0][0] * A[0][1] + A[1][0] * A[1][1] + A[2][0] * A[2][1]
                    ],
                    [
                        A[0][1] * A[0][0] + A[1][1] * A[1][0] + A[2][1] * A[2][0],
                        A[0][1] * A[0][1] + A[1][1] * A[1][1] + A[2][1] * A[2][1]
                    ]
                ];
                const Atb = [
                    A[0][0] * b[0] + A[1][0] * b[1] + A[2][0] * b[2],
                    A[0][1] * b[0] + A[1][1] * b[1] + A[2][1] * b[2]
                ];
                // Solve 2x2 system
                const det = AtA[0][0] * AtA[1][1] - AtA[0][1] * AtA[1][0];
                if (Math.abs(det) < Number.EPSILON) {
                    throw new Error('Singular matrix in Newton iteration');
                }
                const dt1 = (AtA[1][1] * Atb[0] - AtA[0][1] * Atb[1]) / det;
                const dt2 = (-AtA[1][0] * Atb[0] + AtA[0][0] * Atb[1]) / det;
                // Update parameters with damping to ensure convergence
                const damping = 0.5;
                const newT1 = currentT1 + damping * dt1;
                const newT2 = currentT2 + damping * dt2;
                // Check for convergence
                if (Math.abs(newT1 - currentT1) < tolerance && Math.abs(newT2 - currentT2) < tolerance) {
                    return {
                        t1: newT1,
                        t2: newT2,
                        point: curve1.evaluate(newT1),
                        distance: curve1.evaluate(newT1).subtract(curve2.evaluate(newT2)).length()
                    };
                }
                // Update current parameters
                currentT1 = newT1;
                currentT2 = newT2;
                // Ensure parameters stay within domain
                currentT1 = Math.max(curve1.domain.min, Math.min(curve1.domain.max, currentT1));
                currentT2 = Math.max(curve2.domain.min, Math.min(curve2.domain.max, currentT2));
            }
            // If no convergence after max iterations, return null
            return null;
        } catch (error) {
            throw new Error(`Refinement iteration failed: ${error.message}`);
        }
    }
    // Helper method to compute bounding box for a curve
    computeBoundingBox(curve) {
        const numSamples = 20;
        const step = (curve.domain.max - curve.domain.min) / (numSamples - 1);
        let minX = Infinity, minY = Infinity, minZ = Infinity;
        let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
        for (let i = 0; i < numSamples; i++) {
            const t = curve.domain.min + i * step;
            const point = curve.evaluate(t);
            minX = Math.min(minX, point.x);
            minY = Math.min(minY, point.y);
            minZ = Math.min(minZ, point.z);
            maxX = Math.max(maxX, point.x);
            maxY = Math.max(maxY, point.y);
            maxZ = Math.max(maxZ, point.z);
        }
        return {
            min: new Vector3D(minX, minY, minZ),
            max: new Vector3D(maxX, maxY, maxZ)
        };
    }
    // Helper method to compute bounding box for a curve interval
    computeIntervalBoundingBox(curve, t1, t2) {
        const numSamples = 5;
        const step = (t2 - t1) / (numSamples - 1);
        let minX = Infinity, minY = Infinity, minZ = Infinity;
        let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
        for (let i = 0; i < numSamples; i++) {
            const t = t1 + i * step;
            const point = curve.evaluate(t);
            minX = Math.min(minX, point.x);
            minY = Math.min(minY, point.y);
            minZ = Math.min(minZ, point.z);
            maxX = Math.max(maxX, point.x);
            maxY = Math.max(maxY, point.y);
            maxZ = Math.max(maxZ, point.z);
        }
        return {
            min: new Vector3D(minX, minY, minZ),
            max: new Vector3D(maxX, maxY, maxZ)
        };
    }
    // Helper method to check if two bounding boxes intersect
    boundingBoxesIntersect(box1, box2) {
        return !(box1.max.x < box2.min.x || box1.min.x > box2.max.x || box1.max.y < box2.min.y || box1.min.y > box2.max.y || box1.max.z < box2.min.z || box1.min.z > box2.max.z);
    }
    computePatchBoundingBox(surface, patch) {
        const numSamples = 3;
        const points = [];
        // Sample points on the patch
        for (let i = 0; i <= numSamples; i++) {
            const u = patch.u1 + (patch.u2 - patch.u1) * (i / numSamples);
            for (let j = 0; j <= numSamples; j++) {
                const v = patch.v1 + (patch.v2 - patch.v1) * (j / numSamples);
                points.push(surface.evaluate(u, v));
            }
        }
        // Calculate bounding box from sample points
        let minX = Infinity, minY = Infinity, minZ = Infinity;
        let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
        for (const point of points) {
            minX = Math.min(minX, point.x);
            minY = Math.min(minY, point.y);
            minZ = Math.min(minZ, point.z);
            maxX = Math.max(maxX, point.x);
            maxY = Math.max(maxY, point.y);
            maxZ = Math.max(maxZ, point.z);
        }
        return {
            min: new Vector3D(minX, minY, minZ),
            max: new Vector3D(maxX, maxY, maxZ)
        };
    }
    reconstructIntersectionCurves(segments, tolerance) {
        const curves = [];
        const used = new Set();
        for (let i = 0; i < segments.length; i++) {
            if (used.has(i))
                continue;
            const curve = [segments[i]];
            used.add(i);
            // Find connected segments in both directions
            let changed;
            do {
                changed = false;
                // Forward search
                const last = curve[curve.length - 1];
                for (let j = 0; j < segments.length; j++) {
                    if (!used.has(j) && last.point.subtract(segments[j].point).length() < tolerance) {
                        curve.push(segments[j]);
                        used.add(j);
                        changed = true;
                        break;
                    }
                }
                // Backward search
                const first = curve[0];
                for (let j = 0; j < segments.length; j++) {
                    if (!used.has(j) && first.point.subtract(segments[j].point).length() < tolerance) {
                        curve.unshift(segments[j]);
                        used.add(j);
                        changed = true;
                        break;
                    }
                }
            } while (changed);
            curves.push(curve);
        }
        return curves;
    }
}
class IntersectionResult {
    constructor() {
        this.points = [];
        this.curves = [];
        this.parameters = [];
    }
    getIntersectionPoints() {
        // Return empty array if no points are stored
        if (!this.points || !Array.isArray(this.points)) {
            return Object.freeze([]);
        }
        try {
            // Create a deep copy of intersection points to prevent modification
            const intersectionPoints = this.points.map(point => {
                // Verify point is a Vector3D instance
                if (!(point instanceof Vector3D)) {
                    throw new Error('Invalid intersection point data');
                }
                // Create new Vector3D instance
                return new Vector3D(point.x === 0 ? 0 : point.x, // Protect against -0
                    point.y === 0 ? 0 : point.y, point.z === 0 ? 0 : point.z);
            });
            // Validate points
            for (const point of intersectionPoints) {
                if (!Number.isFinite(point.x) || !Number.isFinite(point.y) || !Number.isFinite(point.z)) {
                    throw new Error('Non-finite values in intersection points');
                }
            }
            // Return immutable array of intersection points
            return Object.freeze(intersectionPoints);
        } catch (error) {
            throw new Error(`Failed to get intersection points: ${error.message}`);
        }
    }
    getIntersectionCurves() {
    }
    getParameters() {
    }
    classify() {
    }


    getIntersectionCurves() {
        // Return empty array if no curves are stored
        if (!this.curves || !Array.isArray(this.curves)) {
            return Object.freeze([]);
        }

        try {
            // Create a deep copy of intersection curves
            const intersectionCurves = this.curves.map(curve => {
                // Verify each curve segment is an array
                if (!Array.isArray(curve)) {
                    throw new Error('Invalid curve segment data');
                }

                // Map each point in the curve segment
                const segmentPoints = curve.map(segment => {
                    // Validate segment structure
                    if (!segment || !segment.point || !segment.params) {
                        throw new Error('Invalid intersection curve segment data');
                    }

                    // Verify point is a Vector3D instance
                    if (!(segment.point instanceof Vector3D)) {
                        throw new Error('Invalid point data in curve segment');
                    }

                    // Create deep copy of segment data
                    return {
                        point: new Vector3D(
                            segment.point.x === 0 ? 0 : segment.point.x, // Protect against -0
                            segment.point.y === 0 ? 0 : segment.point.y,
                            segment.point.z === 0 ? 0 : segment.point.z
                        ),
                        params: Object.freeze({
                            surface1: { ...segment.params.surface1 },
                            surface2: { ...segment.params.surface2 }
                        })
                    };
                });

                // Validate all points in segment
                for (const segment of segmentPoints) {
                    if (!Number.isFinite(segment.point.x) ||
                        !Number.isFinite(segment.point.y) ||
                        !Number.isFinite(segment.point.z)) {
                        throw new Error('Non-finite values in intersection curve points');
                    }
                }

                // Return immutable curve segment
                return Object.freeze(segmentPoints);
            });

            // Return immutable array of intersection curves
            return Object.freeze(intersectionCurves);

        } catch (error) {
            throw new Error(`Failed to get intersection curves: ${error.message}`);
        }
    }
}
class GeometryContainer {
    constructor() {
        this.geometries = new Map();
    }
    add(id, geometry) {
        // Input validation
        if (typeof id !== 'string' || id.trim().length === 0) {
            throw new Error('ID must be a non-empty string');
        }
        if (!geometry) {
            throw new Error('Geometry object cannot be null or undefined');
        }
        // Validate geometry type
        if (!(geometry instanceof NurbsCurve) && !(geometry instanceof NurbsSurface) && !(geometry instanceof TrimmedSurface)) {
            throw new Error('Invalid geometry type. Must be NurbsCurve, NurbsSurface, or TrimmedSurface');
        }
        // Check if ID already exists
        if (this.geometries.has(id)) {
            throw new Error(`Geometry with ID "${id}" already exists`);
        }
        try {
            // Create a deep copy of the geometry to prevent external modifications
            let geometryCopy;
            if (geometry instanceof NurbsCurve) {
                geometryCopy = new NurbsCurve([...geometry.controlPoints], new KnotVector([...geometry.knotVector.knots]), geometry.degree);
            } else if (geometry instanceof NurbsSurface) {
                geometryCopy = new NurbsSurface();
                geometryCopy.controlPoints = geometry.controlPoints.map(row => [...row]);
                geometryCopy.knotVectorU = new KnotVector([...geometry.knotVectorU.knots]);
                geometryCopy.knotVectorV = new KnotVector([...geometry.knotVectorV.knots]);
                geometryCopy.degreeU = geometry.degreeU;
                geometryCopy.degreeV = geometry.degreeV;
            } else if (geometry instanceof TrimmedSurface) {
                geometryCopy = new TrimmedSurface();
                geometryCopy.baseSurface = geometry.baseSurface;
                geometryCopy.trimCurves = geometry.trimCurves.map(loop => ({
                    curves: [...loop.curves],
                    isOuter: loop.isOuter
                }));
            }
            // Store the geometry copy
            this.geometries.set(id, {
                type: geometry.constructor.name,
                data: geometryCopy,
                dateAdded: new Date(),
                lastModified: new Date()
            });
            // Return the ID for reference
            return id;
        } catch (error) {
            throw new Error(`Failed to add geometry: ${error.message}`);
        }
    }
    remove(id) {
        // Input validation
        if (typeof id !== 'string' || id.trim().length === 0) {
            throw new Error('ID must be a non-empty string');
        }
        // Check if geometry exists
        if (!this.geometries.has(id)) {
            throw new Error(`Geometry with ID "${id}" does not exist`);
        }
        try {
            // Get the geometry before removal for return value
            const geometry = this.geometries.get(id);
            // Remove the geometry
            const success = this.geometries.delete(id);
            if (!success) {
                throw new Error(`Failed to remove geometry with ID "${id}"`);
            }
            // Return the removed geometry data
            return {
                id: id,
                type: geometry.type,
                data: geometry.data,
                dateAdded: geometry.dateAdded,
                lastModified: geometry.lastModified
            };
        } catch (error) {
            throw new Error(`Error removing geometry: ${error.message}`);
        }
    }
    find(criteria) {
        // Input validation
        if (!criteria || typeof criteria !== 'object') {
            throw new Error('Search criteria must be provided as an object');
        }
        try {
            // Initialize results array
            const results = [];
            // Iterate through all geometries
            for (const [id, geometry] of this.geometries) {
                let matches = true;
                // Check each criterion
                for (const [key, value] of Object.entries(criteria)) {
                    switch (key) {
                        case 'id':
                            matches = matches && id === value;
                            break;
                        case 'type':
                            matches = matches && geometry.type === value;
                            break;
                        case 'dateRange':
                            if (value.from && geometry.dateAdded < value.from) {
                                matches = false;
                            }
                            if (value.to && geometry.dateAdded > value.to) {
                                matches = false;
                            }
                            break;
                        case 'modifiedRange':
                            if (value.from && geometry.lastModified < value.from) {
                                matches = false;
                            }
                            if (value.to && geometry.lastModified > value.to) {
                                matches = false;
                            }
                            break;
                        case 'boundingBox':
                            if (geometry.data instanceof NurbsCurve || geometry.data instanceof NurbsSurface || geometry.data instanceof TrimmedSurface) {
                                const bbox = this.calculateBoundingBox(geometry.data);
                                matches = matches && this.isInBoundingBox(bbox, value);
                            }
                            break;
                        case 'degree':
                            if (geometry.data instanceof NurbsCurve) {
                                matches = matches && geometry.data.degree === value;
                            } else if (geometry.data instanceof NurbsSurface) {
                                matches = matches && ((value.u === undefined || geometry.data.degreeU === value.u) && (value.v === undefined || geometry.data.degreeV === value.v));
                            }
                            break;
                        default:
                            throw new Error(`Invalid search criterion: ${key}`);
                    }
                    // Early exit if no match
                    if (!matches)
                        break;
                }
                // If all criteria match, add to results
                if (matches) {
                    results.push({
                        id: id,
                        type: geometry.type,
                        data: geometry.data,
                        dateAdded: geometry.dateAdded,
                        lastModified: geometry.lastModified
                    });
                }
            }
            return Object.freeze(results);
        } catch (error) {
            throw new Error(`Search failed: ${error.message}`);
        }
    }
    transform(id, transformationMatrix) {
        // Input validation
        if (typeof id !== 'string' || id.trim().length === 0) {
            throw new Error('ID must be a non-empty string');
        }
        if (!(transformationMatrix instanceof Matrix4x4)) {
            throw new Error('Transformation matrix must be an instance of Matrix4x4');
        }
        // Check if geometry exists
        if (!this.geometries.has(id)) {
            throw new Error(`Geometry with ID "${id}" does not exist`);
        }
        const geometry = this.geometries.get(id);
        let transformedGeometry;
        try {
            if (geometry.type === 'NurbsCurve') {
                // Transform curve control points
                const transformedControlPoints = geometry.data.controlPoints.map(cp => {
                    const [x, y, z] = cp.position();
                    const w = cp.weight();
                    // Create homogeneous point
                    const point = [
                        x * w,
                        y * w,
                        z * w,
                        w
                    ];
                    // Apply transformation
                    const transformed = transformationMatrix.multiply(new Matrix4x4([
                        point[0],
                        point[1],
                        point[2],
                        point[3],
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0
                    ]));
                    // Extract transformed coordinates
                    const transformedPoint = transformed.elements;
                    const newW = transformedPoint[3];
                    if (Math.abs(newW) < Number.EPSILON) {
                        throw new Error('Invalid transformation: results in zero weight');
                    }
                    // Create new control point
                    return new ControlPoint(transformedPoint[0] / newW, transformedPoint[1] / newW, transformedPoint[2] / newW, newW);
                });
                // Create new curve with transformed control points
                transformedGeometry = new NurbsCurve(transformedControlPoints, geometry.data.knotVector, geometry.data.degree);
            } else if (geometry.type === 'NurbsSurface') {
                // Transform surface control points
                const transformedControlPoints = geometry.data.controlPoints.map(row => row.map(cp => {
                    const [x, y, z] = cp.position();
                    const w = cp.weight();
                    // Create homogeneous point
                    const point = [
                        x * w,
                        y * w,
                        z * w,
                        w
                    ];
                    // Apply transformation
                    const transformed = transformationMatrix.multiply(new Matrix4x4([
                        point[0],
                        point[1],
                        point[2],
                        point[3],
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0
                    ]));
                    // Extract transformed coordinates
                    const transformedPoint = transformed.elements;
                    const newW = transformedPoint[3];
                    if (Math.abs(newW) < Number.EPSILON) {
                        throw new Error('Invalid transformation: results in zero weight');
                    }
                    // Create new control point
                    return new ControlPoint(transformedPoint[0] / newW, transformedPoint[1] / newW, transformedPoint[2] / newW, newW);
                }));
                // Create new surface with transformed control points
                transformedGeometry = new NurbsSurface();
                transformedGeometry.controlPoints = transformedControlPoints;
                transformedGeometry.knotVectorU = geometry.data.knotVectorU;
                transformedGeometry.knotVectorV = geometry.data.knotVectorV;
                transformedGeometry.degreeU = geometry.data.degreeU;
                transformedGeometry.degreeV = geometry.data.degreeV;
            } else if (geometry.type === 'TrimmedSurface') {
                // Transform base surface
                const transformedBaseSurface = new NurbsSurface();
                transformedBaseSurface.controlPoints = geometry.data.baseSurface.controlPoints.map(row => row.map(cp => {
                    const [x, y, z] = cp.position();
                    const w = cp.weight();
                    const point = [
                        x * w,
                        y * w,
                        z * w,
                        w
                    ];
                    const transformed = transformationMatrix.multiply(new Matrix4x4([
                        point[0],
                        point[1],
                        point[2],
                        point[3],
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0
                    ]));
                    const transformedPoint = transformed.elements;
                    const newW = transformedPoint[3];
                    if (Math.abs(newW) < Number.EPSILON) {
                        throw new Error('Invalid transformation: results in zero weight');
                    }
                    return new ControlPoint(transformedPoint[0] / newW, transformedPoint[1] / newW, transformedPoint[2] / newW, newW);
                }));
                transformedBaseSurface.knotVectorU = geometry.data.baseSurface.knotVectorU;
                transformedBaseSurface.knotVectorV = geometry.data.baseSurface.knotVectorV;
                transformedBaseSurface.degreeU = geometry.data.baseSurface.degreeU;
                transformedBaseSurface.degreeV = geometry.data.baseSurface.degreeV;
                // Create new trimmed surface
                transformedGeometry = new TrimmedSurface();
                transformedGeometry.baseSurface = transformedBaseSurface;
                transformedGeometry.trimCurves = geometry.data.trimCurves;
            }
            // Update geometry in container
            this.geometries.set(id, {
                type: geometry.type,
                data: transformedGeometry,
                dateAdded: geometry.dateAdded,
                lastModified: new Date()
            });
            return transformedGeometry;
        } catch (error) {
            throw new Error(`Transformation failed: ${error.message}`);
        }
    }
    validate() {
        const validationResults = {
            isValid: true,
            errors: [],
            warnings: []
        };
        try {
            // Validate each geometry in the container
            for (const [id, geometry] of this.geometries) {
                // Check ID format
                if (typeof id !== 'string' || id.trim().length === 0) {
                    validationResults.errors.push(`Invalid ID format for geometry: ${id}`);
                    validationResults.isValid = false;
                    continue;
                }
                // Check geometry object structure
                if (!geometry || !geometry.type || !geometry.data) {
                    validationResults.errors.push(`Invalid geometry structure for ID: ${id}`);
                    validationResults.isValid = false;
                    continue;
                }
                // Check geometry type
                if (![
                    'NurbsCurve',
                    'NurbsSurface',
                    'TrimmedSurface'
                ].includes(geometry.type)) {
                    validationResults.errors.push(`Invalid geometry type for ID: ${id}`);
                    validationResults.isValid = false;
                    continue;
                }
                // Type-specific validation
                switch (geometry.type) {
                    case 'NurbsCurve':
                        this.validateNurbsCurve(id, geometry.data, validationResults);
                        break;
                    case 'NurbsSurface':
                        this.validateNurbsSurface(id, geometry.data, validationResults);
                        break;
                    case 'TrimmedSurface':
                        this.validateTrimmedSurface(id, geometry.data, validationResults);
                        break;
                }
                // Validate timestamps
                if (!(geometry.dateAdded instanceof Date) || isNaN(geometry.dateAdded.getTime())) {
                    validationResults.errors.push(`Invalid dateAdded for ID: ${id}`);
                    validationResults.isValid = false;
                }
                if (!(geometry.lastModified instanceof Date) || isNaN(geometry.lastModified.getTime())) {
                    validationResults.errors.push(`Invalid lastModified for ID: ${id}`);
                    validationResults.isValid = false;
                }
                // Check for future dates
                const now = new Date();
                if (geometry.dateAdded > now || geometry.lastModified > now) {
                    validationResults.warnings.push(`Future timestamps detected for ID: ${id}`);
                }
            }
            // Check for duplicate IDs
            const ids = Array.from(this.geometries.keys());
            const uniqueIds = new Set(ids);
            if (ids.length !== uniqueIds.size) {
                validationResults.errors.push('Duplicate geometry IDs detected');
                validationResults.isValid = false;
            }
            return Object.freeze(validationResults);
        } catch (error) {
            validationResults.errors.push(`Validation failed: ${error.message}`);
            validationResults.isValid = false;
            return Object.freeze(validationResults);
        }
    }
    // Helper method to calculate bounding box
    calculateBoundingBox(geometry) {
        const bbox = {
            min: new Vector3D(Infinity, Infinity, Infinity),
            max: new Vector3D(-Infinity, -Infinity, -Infinity)
        };
        if (geometry instanceof NurbsCurve) {
            // Sample curve at regular intervals
            const numSamples = 50;
            const step = (geometry.domain.max - geometry.domain.min) / (numSamples - 1);
            for (let i = 0; i < numSamples; i++) {
                const point = geometry.evaluate(geometry.domain.min + i * step);
                this.updateBoundingBox(bbox, point);
            }
        } else if (geometry instanceof NurbsSurface || geometry instanceof TrimmedSurface) {
            const surface = geometry instanceof TrimmedSurface ? geometry.baseSurface : geometry;
            // Sample surface at regular intervals
            const numSamplesU = 20;
            const numSamplesV = 20;
            const stepU = (surface.knotVectorU.domain.max - surface.knotVectorU.domain.min) / (numSamplesU - 1);
            const stepV = (surface.knotVectorV.domain.max - surface.knotVectorV.domain.min) / (numSamplesV - 1);
            for (let i = 0; i < numSamplesU; i++) {
                for (let j = 0; j < numSamplesV; j++) {
                    const u = surface.knotVectorU.domain.min + i * stepU;
                    const v = surface.knotVectorV.domain.min + j * stepV;
                    const point = surface.evaluate(u, v);
                    this.updateBoundingBox(bbox, point);
                }
            }
        }
        return bbox;
    }
    // Helper method to update bounding box with new point
    updateBoundingBox(bbox, point) {
        bbox.min.x = Math.min(bbox.min.x, point.x);
        bbox.min.y = Math.min(bbox.min.y, point.y);
        bbox.min.z = Math.min(bbox.min.z, point.z);
        bbox.max.x = Math.max(bbox.max.x, point.x);
        bbox.max.y = Math.max(bbox.max.y, point.y);
        bbox.max.z = Math.max(bbox.max.z, point.z);
    }
    // Helper method to check if bounding box is within search bounds
    isInBoundingBox(bbox, searchBox) {
        return !(bbox.min.x > searchBox.max.x || bbox.max.x < searchBox.min.x || bbox.min.y > searchBox.max.y || bbox.max.y < searchBox.min.y || bbox.min.z > searchBox.max.z || bbox.max.z < searchBox.min.z);
    }
    validateNurbsCurve(id, curve, results) {
        if (!(curve instanceof NurbsCurve)) {
            results.errors.push(`Invalid NURBS curve instance for ID: ${id}`);
            results.isValid = false;
            return;
        }
        // Validate control points
        if (!Array.isArray(curve.controlPoints) || curve.controlPoints.length < curve.degree + 1) {
            results.errors.push(`Invalid control points for curve ID: ${id}`);
            results.isValid = false;
        }
        // Validate knot vector
        if (!(curve.knotVector instanceof KnotVector)) {
            results.errors.push(`Invalid knot vector for curve ID: ${id}`);
            results.isValid = false;
        }
        // Validate degree
        if (!Number.isInteger(curve.degree) || curve.degree < 1) {
            results.errors.push(`Invalid degree for curve ID: ${id}`);
            results.isValid = false;
        }
        // Check for degenerate curves
        try {
            const start = curve.evaluate(curve.domain.min);
            const end = curve.evaluate(curve.domain.max);
            if (start.subtract(end).length() < Number.EPSILON) {
                results.warnings.push(`Potentially degenerate curve detected for ID: ${id}`);
            }
        } catch (error) {
            results.errors.push(`Curve evaluation failed for ID: ${id}: ${error.message}`);
            results.isValid = false;
        }
    }
    validateNurbsSurface(id, surface, results) {
        if (!(surface instanceof NurbsSurface)) {
            results.errors.push(`Invalid NURBS surface instance for ID: ${id}`);
            results.isValid = false;
            return;
        }
        // Validate control point network
        if (!Array.isArray(surface.controlPoints) || surface.controlPoints.length === 0) {
            results.errors.push(`Invalid control point network for surface ID: ${id}`);
            results.isValid = false;
        }
        // Check control point network regularity
        const width = surface.controlPoints[0].length;
        if (!surface.controlPoints.every(row => Array.isArray(row) && row.length === width)) {
            results.errors.push(`Irregular control point network for surface ID: ${id}`);
            results.isValid = false;
        }
        // Validate knot vectors
        if (!(surface.knotVectorU instanceof KnotVector) || !(surface.knotVectorV instanceof KnotVector)) {
            results.errors.push(`Invalid knot vectors for surface ID: ${id}`);
            results.isValid = false;
        }
        // Validate degrees
        if (!Number.isInteger(surface.degreeU) || !Number.isInteger(surface.degreeV) || surface.degreeU < 1 || surface.degreeV < 1) {
            results.errors.push(`Invalid degrees for surface ID: ${id}`);
            results.isValid = false;
        }
        // Check for degenerate surfaces
        try {
            const area = this.calculateApproximateArea(surface);
            if (area < Number.EPSILON) {
                results.warnings.push(`Potentially degenerate surface detected for ID: ${id}`);
            }
        } catch (error) {
            results.errors.push(`Surface evaluation failed for ID: ${id}: ${error.message}`);
            results.isValid = false;
        }
    }
    validateTrimmedSurface(id, surface, results) {
        if (!(surface instanceof TrimmedSurface)) {
            results.errors.push(`Invalid trimmed surface instance for ID: ${id}`);
            results.isValid = false;
            return;
        }
        // Validate base surface
        if (!(surface.baseSurface instanceof NurbsSurface)) {
            results.errors.push(`Invalid base surface for trimmed surface ID: ${id}`);
            results.isValid = false;
        }
        // Validate trim curves
        if (!Array.isArray(surface.trimCurves)) {
            results.errors.push(`Invalid trim curves for surface ID: ${id}`);
            results.isValid = false;
            return;
        }
        // Check trim curve validity
        for (let i = 0; i < surface.trimCurves.length; i++) {
            const trimLoop = surface.trimCurves[i];
            if (!Array.isArray(trimLoop.curves) || !trimLoop.hasOwnProperty('isOuter')) {
                results.errors.push(`Invalid trim loop structure at index ${i} for surface ID: ${id}`);
                results.isValid = false;
                continue;
            }
            // Verify trim curves form closed loops
            const curves = trimLoop.curves;
            for (let j = 0; j < curves.length; j++) {
                if (!(curves[j] instanceof NurbsCurve)) {
                    results.errors.push(`Invalid trim curve at index ${j} in loop ${i} for surface ID: ${id}`);
                    results.isValid = false;
                }
            }
        }
        // Check for valid outer boundary
        const hasOuterBoundary = surface.trimCurves.some(loop => loop.isOuter);
        if (!hasOuterBoundary && surface.trimCurves.length > 0) {
            results.errors.push(`Missing outer boundary for trimmed surface ID: ${id}`);
            results.isValid = false;
        }
    }
    calculateApproximateArea(surface) {
        const numSamples = 10;
        let area = 0;
        const uMin = surface.knotVectorU.knots[surface.degreeU];
        const uMax = surface.knotVectorU.knots[surface.knotVectorU.length - surface.degreeU - 1];
        const vMin = surface.knotVectorV.knots[surface.degreeV];
        const vMax = surface.knotVectorV.knots[surface.knotVectorV.length - surface.degreeV - 1];
        const du = (uMax - uMin) / numSamples;
        const dv = (vMax - vMin) / numSamples;
        for (let i = 0; i < numSamples; i++) {
            for (let j = 0; j < numSamples; j++) {
                const u = uMin + i * du;
                const v = vMin + j * dv;
                const p00 = surface.evaluate(u, v);
                const p10 = surface.evaluate(u + du, v);
                const p01 = surface.evaluate(u, v + dv);
                const v1 = p10.subtract(p00);
                const v2 = p01.subtract(p00);
                area += v1.cross(v2).length() / 2;
            }
        }
        return area;
    }
}
class FileIO {
    constructor() {
    }
    readIGES() {
    }
    writeIGES() {
    }
    readSTEP() {
    }
    writeSTEP() {
    }
}
class Tessellator {
    constructor() {
        this.tolerance = defaultTolerance;
    }
    tessellateNurbsCurve(curve, tolerance = this.tolerance) {
        // Input validation
        if (!(curve instanceof NurbsCurve)) {
            throw new Error('Input must be a NURBS curve');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        // Initialize points array with start point
        const points = [curve.evaluate(curve.domain.min)];
        const parameters = [curve.domain.min];
        // Recursive subdivision function
        const subdivide = (u1, u2) => {
            // Evaluate curve at endpoints and midpoint
            const p1 = curve.evaluate(u1);
            const p2 = curve.evaluate(u2);
            const umid = (u1 + u2) / 2;
            const pmid = curve.evaluate(umid);
            // Calculate chord length
            const chordLength = p1.subtract(p2).length();
            // Calculate maximum deviation from chord
            const deviation = pmid.subtract(p1.add(p2.subtract(p1).multiply(0.5))).length();
            // Calculate local curvature at midpoint
            let curvature = 0;
            try {
                const d1 = curve.derivative(umid, 1);
                const d2 = curve.derivative(umid, 2);
                const numerator = d1.cross(d2).length();
                const denominator = Math.pow(d1.length(), 3);
                if (denominator > Number.EPSILON) {
                    curvature = numerator / denominator;
                }
            } catch (error) {
                // If curvature calculation fails, use deviation only
                curvature = 0;
            }
            // Determine if subdivision is needed based on tolerance and curvature
            const maxDeviation = Math.max(deviation, curvature * chordLength * chordLength / 8);
            if (maxDeviation > tolerance && Math.abs(u2 - u1) > Number.EPSILON) {
                // Subdivide further
                subdivide(u1, umid);
                // Add midpoint
                points.push(pmid);
                parameters.push(umid);
                subdivide(umid, u2);
            }
        };
        try {
            // Perform adaptive subdivision
            subdivide(curve.domain.min, curve.domain.max);
            // Add end point
            points.push(curve.evaluate(curve.domain.max));
            parameters.push(curve.domain.max);
            // Remove any duplicate points within tolerance
            const filteredPoints = [points[0]];
            const filteredParameters = [parameters[0]];
            for (let i = 1; i < points.length; i++) {
                if (points[i].subtract(filteredPoints[filteredPoints.length - 1]).length() > tolerance) {
                    filteredPoints.push(points[i]);
                    filteredParameters.push(parameters[i]);
                }
            }
            // Ensure the end point is included
            const lastPoint = curve.evaluate(curve.domain.max);
            if (lastPoint.subtract(filteredPoints[filteredPoints.length - 1]).length() > tolerance) {
                filteredPoints.push(lastPoint);
                filteredParameters.push(curve.domain.max);
            }
            // Return tessellation result
            return {
                points: Object.freeze(filteredPoints),
                parameters: Object.freeze(filteredParameters)
            };
        } catch (error) {
            throw new Error(`Tessellation failed: ${error.message}`);
        }
    }
    tessellateNurbsSurface(surface, tolerance = this.tolerance) {
        // Input validation
        if (!(surface instanceof NurbsSurface)) {
            throw new Error('Input must be a NURBS surface');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        // Get surface domain
        const uDomain = {
            min: surface.knotVectorU.knots[surface.degreeU],
            max: surface.knotVectorU.knots[surface.knotVectorU.length - surface.degreeU - 1]
        };
        const vDomain = {
            min: surface.knotVectorV.knots[surface.degreeV],
            max: surface.knotVectorV.knots[surface.knotVectorV.length - surface.degreeV - 1]
        };
        // Initialize points and parameters arrays
        const points = [];
        const parameters = [];
        const triangles = [];
        // Recursive subdivision function for a surface patch
        const subdividePatch = (u1, u2, v1, v2, depth = 0) => {
            const maxDepth = 10;
            // Maximum recursion depth to prevent infinite recursion
            if (depth > maxDepth) {
                return;
            }
            // Evaluate surface at corners and mid-points
            const p00 = surface.evaluate(u1, v1);
            const p10 = surface.evaluate(u2, v1);
            const p01 = surface.evaluate(u1, v2);
            const p11 = surface.evaluate(u2, v2);
            const uMid = (u1 + u2) / 2;
            const vMid = (v1 + v2) / 2;
            const pMidMid = surface.evaluate(uMid, vMid);
            // Calculate maximum deviation from bilinear interpolation
            const interpolated = p00.add(p10).add(p01).add(p11).multiply(0.25);
            const deviation = pMidMid.subtract(interpolated).length();
            // Calculate local curvature at midpoint
            let maxCurvature = 0;
            try {
                const curvatureInfo = surface.curvature(uMid, vMid);
                maxCurvature = Math.max(Math.abs(curvatureInfo.principalCurvatures.k1), Math.abs(curvatureInfo.principalCurvatures.k2));
            } catch (error) {
                // If curvature calculation fails, use deviation only
                maxCurvature = 0;
            }
            // Calculate patch size
            const patchSize = Math.max(p00.subtract(p10).length(), p00.subtract(p01).length());
            // Determine if subdivision is needed based on tolerance and curvature
            const maxDeviation = Math.max(deviation, maxCurvature * patchSize * patchSize / 8);
            if (maxDeviation > tolerance && Math.abs(u2 - u1) > Number.EPSILON && Math.abs(v2 - v1) > Number.EPSILON) {
                // Subdivide patch into four smaller patches
                subdividePatch(u1, uMid, v1, vMid, depth + 1);
                subdividePatch(uMid, u2, v1, vMid, depth + 1);
                subdividePatch(u1, uMid, vMid, v2, depth + 1);
                subdividePatch(uMid, u2, vMid, v2, depth + 1);
            } else {
                // Add points and create triangles for this patch
                const index = points.length;
                points.push(p00, p10, p01, p11);
                parameters.push({
                    u: u1,
                    v: v1
                }, {
                    u: u2,
                    v: v1
                }, {
                    u: u1,
                    v: v2
                }, {
                    u: u2,
                    v: v2
                });
                // Create two triangles for the patch
                triangles.push([
                    index,
                    index + 1,
                    index + 2
                ], [
                    index + 1,
                    index + 3,
                    index + 2
                ]);
            }
        };
        try {
            // Start recursive subdivision with entire surface
            subdividePatch(uDomain.min, uDomain.max, vDomain.min, vDomain.max);
            // Remove duplicate points within tolerance
            const uniquePoints = [];
            const uniqueParameters = [];
            const indexMap = new Map();
            for (let i = 0; i < points.length; i++) {
                let foundMatch = false;
                for (let j = 0; j < uniquePoints.length; j++) {
                    if (points[i].subtract(uniquePoints[j]).length() < tolerance) {
                        indexMap.set(i, j);
                        foundMatch = true;
                        break;
                    }
                }
                if (!foundMatch) {
                    indexMap.set(i, uniquePoints.length);
                    uniquePoints.push(points[i]);
                    uniqueParameters.push(parameters[i]);
                }
            }
            // Remap triangle indices
            const remappedTriangles = triangles.map(triangle => triangle.map(index => indexMap.get(index)));
            // Return tessellation result
            return {
                points: Object.freeze(uniquePoints),
                parameters: Object.freeze(uniqueParameters),
                triangles: Object.freeze(remappedTriangles)
            };
        } catch (error) {
            throw new Error(`Surface tessellation failed: ${error.message}`);
        }
    }
    tessellateTrimmedSurface(surface, tolerance = this.tolerance) {
        // Input validation
        if (!(surface instanceof TrimmedSurface)) {
            throw new Error('Input must be a TrimmedSurface instance');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        try {
            // First, tessellate the base surface
            const baseSurfaceTess = this.tessellateNurbsSurface(surface.baseSurface, tolerance);
            if (!baseSurfaceTess || !baseSurfaceTess.points || !baseSurfaceTess.parameters) {
                throw new Error('Base surface tessellation failed');
            }
            // Create arrays to store valid points and triangles
            const validPoints = [];
            const validParameters = [];
            const validTriangles = [];
            // Map to store original point indices to new indices
            const indexMap = new Map();
            // Check each point against trim curves
            for (let i = 0; i < baseSurfaceTess.points.length; i++) {
                const param = baseSurfaceTess.parameters[i];
                if (surface.isPointInside(param.u, param.v)) {
                    indexMap.set(i, validPoints.length);
                    validPoints.push(baseSurfaceTess.points[i]);
                    validParameters.push(param);
                }
            }
            // Process triangles, keeping only those with all points inside
            for (const triangle of baseSurfaceTess.triangles) {
                const newTriangle = triangle.map(index => indexMap.get(index));
                if (newTriangle.every(index => index !== undefined)) {
                    validTriangles.push(newTriangle);
                }
            }
            // Tessellate trim curves
            const tessellatedTrimCurves = [];
            for (const trimLoop of surface.trimCurves) {
                const loopCurves = [];
                for (const curve of trimLoop.curves) {
                    const curveTess = this.tessellateNurbsCurve(curve, tolerance);
                    loopCurves.push({
                        points: curveTess.points,
                        parameters: curveTess.parameters,
                        isOuter: trimLoop.isOuter
                    });
                }
                tessellatedTrimCurves.push(loopCurves);
            }
            // Add points along trim curves to ensure proper boundary representation
            const boundaryPoints = new Set();
            const boundaryParameters = new Set();
            for (const trimLoop of tessellatedTrimCurves) {
                for (const curve of trimLoop) {
                    for (let i = 0; i < curve.points.length; i++) {
                        const point = curve.points[i];
                        const param = curve.parameters[i];
                        // Check if point is already included (within tolerance)
                        let foundMatch = false;
                        for (let j = 0; j < validPoints.length; j++) {
                            if (point.subtract(validPoints[j]).length() < tolerance) {
                                foundMatch = true;
                                break;
                            }
                        }
                        if (!foundMatch) {
                            boundaryPoints.add(point);
                            boundaryParameters.add({
                                u: param.u,
                                v: param.v
                            });
                        }
                    }
                }
            }
            // Add boundary points to valid points array
            const boundaryPointsArray = Array.from(boundaryPoints);
            const boundaryParamsArray = Array.from(boundaryParameters);
            for (let i = 0; i < boundaryPointsArray.length; i++) {
                validPoints.push(boundaryPointsArray[i]);
                validParameters.push(boundaryParamsArray[i]);
            }
            // Return tessellation result
            return {
                points: Object.freeze(validPoints),
                parameters: Object.freeze(validParameters),
                triangles: Object.freeze(validTriangles),
                trimCurves: Object.freeze(tessellatedTrimCurves),
                tolerance: tolerance
            };
        } catch (error) {
            throw new Error(`Tessellation of trimmed surface failed: ${error.message}`);
        }
    }
    generateKnotLines(surface, tolerance = this.tolerance) {
        // Input validation
        if (!(surface instanceof NurbsSurface)) {
            throw new Error('Input must be a NURBS surface');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        const knotLinesU = [];
        const knotLinesV = [];
        try {
            // Get unique knots in U direction
            const uniqueKnotsU = Array.from(new Set(surface.knotVectorU.knots)).filter(knot => {
                const min = surface.knotVectorU.knots[surface.degreeU];
                const max = surface.knotVectorU.knots[surface.knotVectorU.length - surface.degreeU - 1];
                return knot >= min && knot <= max;
            });
            // Get unique knots in V direction
            const uniqueKnotsV = Array.from(new Set(surface.knotVectorV.knots)).filter(knot => {
                const min = surface.knotVectorV.knots[surface.degreeV];
                const max = surface.knotVectorV.knots[surface.knotVectorV.length - surface.degreeV - 1];
                return knot >= min && knot <= max;
            });
            // Generate U-direction knot lines
            for (const u of uniqueKnotsU) {
                const knotLine = [];
                const vMin = surface.knotVectorV.knots[surface.degreeV];
                const vMax = surface.knotVectorV.knots[surface.knotVectorV.length - surface.degreeV - 1];
                // Initial points at ends of V domain
                const startPoint = surface.evaluate(u, vMin);
                const endPoint = surface.evaluate(u, vMax);
                knotLine.push(startPoint);
                // Adaptive sampling based on curvature
                const subdivide = (v1, v2) => {
                    const vMid = (v1 + v2) / 2;
                    const p1 = surface.evaluate(u, v1);
                    const p2 = surface.evaluate(u, v2);
                    const pMid = surface.evaluate(u, vMid);
                    // Calculate deviation from linear interpolation
                    const interpolated = p1.add(p2.subtract(p1).multiply(0.5));
                    const deviation = pMid.subtract(interpolated).length();
                    // Calculate local curvature
                    let curvature = 0;
                    try {
                        const du = surface.derivative(u, vMid, 0, 1);
                        const d2u = surface.derivative(u, vMid, 0, 2);
                        const numerator = du.cross(d2u).length();
                        const denominator = Math.pow(du.length(), 3);
                        if (denominator > Number.EPSILON) {
                            curvature = numerator / denominator;
                        }
                    } catch (error) {
                        // If curvature calculation fails, use deviation only
                        curvature = 0;
                    }
                    // Determine if subdivision is needed
                    const chordLength = p1.subtract(p2).length();
                    const maxDeviation = Math.max(deviation, curvature * chordLength * chordLength / 8);
                    if (maxDeviation > tolerance && Math.abs(v2 - v1) > Number.EPSILON) {
                        subdivide(v1, vMid);
                        knotLine.push(pMid);
                        subdivide(vMid, v2);
                    }
                };
                // Perform adaptive subdivision
                subdivide(vMin, vMax);
                knotLine.push(endPoint);
                knotLinesU.push(Object.freeze(knotLine));
            }
            // Generate V-direction knot lines
            for (const v of uniqueKnotsV) {
                const knotLine = [];
                const uMin = surface.knotVectorU.knots[surface.degreeU];
                const uMax = surface.knotVectorU.knots[surface.knotVectorU.length - surface.degreeU - 1];
                // Initial points at ends of U domain
                const startPoint = surface.evaluate(uMin, v);
                const endPoint = surface.evaluate(uMax, v);
                knotLine.push(startPoint);
                // Adaptive sampling based on curvature
                const subdivide = (u1, u2) => {
                    const uMid = (u1 + u2) / 2;
                    const p1 = surface.evaluate(u1, v);
                    const p2 = surface.evaluate(u2, v);
                    const pMid = surface.evaluate(uMid, v);
                    // Calculate deviation from linear interpolation
                    const interpolated = p1.add(p2.subtract(p1).multiply(0.5));
                    const deviation = pMid.subtract(interpolated).length();
                    // Calculate local curvature
                    let curvature = 0;
                    try {
                        const dv = surface.derivative(uMid, v, 1, 0);
                        const d2v = surface.derivative(uMid, v, 2, 0);
                        const numerator = dv.cross(d2v).length();
                        const denominator = Math.pow(dv.length(), 3);
                        if (denominator > Number.EPSILON) {
                            curvature = numerator / denominator;
                        }
                    } catch (error) {
                        // If curvature calculation fails, use deviation only
                        curvature = 0;
                    }
                    // Determine if subdivision is needed
                    const chordLength = p1.subtract(p2).length();
                    const maxDeviation = Math.max(deviation, curvature * chordLength * chordLength / 8);
                    if (maxDeviation > tolerance && Math.abs(u2 - u1) > Number.EPSILON) {
                        subdivide(u1, uMid);
                        knotLine.push(pMid);
                        subdivide(uMid, u2);
                    }
                };
                // Perform adaptive subdivision
                subdivide(uMin, uMax);
                knotLine.push(endPoint);
                knotLinesV.push(Object.freeze(knotLine));
            }
            // Return both sets of knot lines
            return {
                uLines: Object.freeze(knotLinesU),
                vLines: Object.freeze(knotLinesV)
            };
        } catch (error) {
            throw new Error(`Failed to generate knot lines: ${error.message}`);
        }
    }
}


class IntersectionResult {
    // ... other methods ...

    getIntersectionCurves() {
        // Return empty array if no curves are stored
        if (!this.curves || !Array.isArray(this.curves)) {
            return Object.freeze([]);
        }

        try {
            // Create a deep copy of intersection curves
            const intersectionCurves = this.curves.map(curve => {
                // Verify each curve segment is an array
                if (!Array.isArray(curve)) {
                    throw new Error('Invalid curve segment data');
                }

                // Map each point in the curve segment
                const segmentPoints = curve.map(segment => {
                    // Validate segment structure
                    if (!segment || !segment.point || !segment.params) {
                        throw new Error('Invalid intersection curve segment data');
                    }

                    // Verify point is a Vector3D instance
                    if (!(segment.point instanceof Vector3D)) {
                        throw new Error('Invalid point data in curve segment');
                    }

                    // Create deep copy of segment data
                    return {
                        point: new Vector3D(
                            segment.point.x === 0 ? 0 : segment.point.x, // Protect against -0
                            segment.point.y === 0 ? 0 : segment.point.y,
                            segment.point.z === 0 ? 0 : segment.point.z
                        ),
                        params: Object.freeze({
                            surface1: { ...segment.params.surface1 },
                            surface2: { ...segment.params.surface2 }
                        })
                    };
                });

                // Validate all points in segment
                for (const segment of segmentPoints) {
                    if (!Number.isFinite(segment.point.x) ||
                        !Number.isFinite(segment.point.y) ||
                        !Number.isFinite(segment.point.z)) {
                        throw new Error('Non-finite values in intersection curve points');
                    }
                }

                // Return immutable curve segment
                return Object.freeze(segmentPoints);
            });

            // Return immutable array of intersection curves
            return Object.freeze(intersectionCurves);

        } catch (error) {
            throw new Error(`Failed to get intersection curves: ${error.message}`);
        }
    }
}








