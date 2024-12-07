const defaultTolerance = 0.000001;
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
            throw new Error(`Failed to calculate derivative: ${ error.message }`);
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
    clone() {
        try {
            // Create new instance
            const clonedBasis = new BasisFunction();
            // Copy degree
            clonedBasis.degree = this.degree;
            // Clone knot vector if it exists
            if (this.knotVector) {
                // KnotVector class should have its own clone method
                clonedBasis.knotVector = this.knotVector.clone();
            }
            return clonedBasis;
        } catch (error) {
            throw new Error(`Failed to clone BasisFunction: ${ error.message }`);
        }
    }
}
class BrepBuilder {
    constructor() {
        // Initialize any required properties for the BrepBuilder
        this.shells = [];
    }
    createFaceFromSurface(surface) {
        // Input validation
        if (!(surface instanceof NurbsSurface)) {
            throw new Error('Parameter must be a NurbsSurface instance');
        }
        try {
            // Create a new face from the provided surface
            const face = new TopologyFace();
            face.surface = surface;
            // Extract edges from the surface's control points to define the face boundaries
            const edges = [];
            const uCount = surface.controlPoints.length;
            const vCount = surface.controlPoints[0].length;
            // Create edges from the control points
            for (let i = 0; i < uCount; i++) {
                for (let j = 0; j < vCount - 1; j++) {
                    const curveEdge = new TopologyEdge();
                    curveEdge.curve = new NurbsCurve([
                        surface.controlPoints[i][j],
                        surface.controlPoints[i][j + 1]
                    ], surface.knotVectorU, surface.degreeU);
                    edges.push(curveEdge);
                    face.edges.push(curveEdge);
                }
            }
            for (let j = 0; j < vCount; j++) {
                for (let i = 0; i < uCount - 1; i++) {
                    const curveEdge = new TopologyEdge();
                    curveEdge.curve = new NurbsCurve([
                        surface.controlPoints[i][j],
                        surface.controlPoints[i + 1][j]
                    ], surface.knotVectorV, surface.degreeV);
                    edges.push(curveEdge);
                    face.edges.push(curveEdge);
                }
            }
            // Validate face with edges
            face.edges = Array.from(new Set(edges));
            // Ensure unique edges
            face.validateBounds();
            // Validate the created face
            return face;
        } // Return the newly created face
        catch (error) {
            throw new Error(`Failed to create face from surface: ${ error.message }`);
        }
    }
    createEdgeFromCurve(curve) {
        // Input validation
        if (!(curve instanceof NurbsCurve)) {
            throw new Error('Parameter must be a NurbsCurve instance');
        }
        try {
            // Create a new edge from the given curve
            const edge = new TopologyEdge();
            edge.curve = curve.clone();
            // Clone to ensure uniqueness
            edge.startVertex = new TopologyVertex({ position: curve.evaluate(curve.domain.min) });
            edge.endVertex = new TopologyVertex({ position: curve.evaluate(curve.domain.max) });
            // Update vertex edge lists
            edge.startVertex.edges.push(edge);
            edge.endVertex.edges.push(edge);
            // Return the newly created edge
            return edge;
        } catch (error) {
            throw new Error(`Failed to create edge from curve: ${ error.message }`);
        }
    }
    createVertex(position) {
        // Input validation
        if (!(position instanceof Vector3D)) {
            throw new Error('Position must be a Vector3D instance');
        }
        if (!Number.isFinite(position.x) || !Number.isFinite(position.y) || !Number.isFinite(position.z)) {
            throw new Error('Position coordinates must be finite numbers');
        }
        try {
            // Create a new TopologyVertex instance at the specified position
            const newVertex = new TopologyVertex({ position: position.clone() });
            // Add to the internal vertices collection in the BrepBuilder
            this.vertices.push(newVertex);
            return newVertex;
        } catch (error) {
            throw new Error(`Failed to create vertex: ${ error.message }`);
        }
    }
    createShell() {
        try {
            const newShell = new Shell();
            // Add the new shell to the internal collection
            this.shells.push(newShell);
            return newShell;
        } catch (error) {
            throw new Error(`Failed to create shell: ${ error.message }`);
        }
    }
    createSolid() {
        // Create a new Solid instance
        const newSolid = new Solid();
        // Add the current shell to the new solid
        this.shells.forEach(shell => {
            newSolid.shells.push(shell);
        });
        // Return the newly created solid
        return newSolid;
    }
    finalize() {
        try {
            if (this.shells.length === 0) {
                throw new Error('No shells available to finalize the BREP structure');
            }
            // Create a new Solid instance to represent the finalized BREP
            const finalizedSolid = new Solid();
            // Add all shells to the finalized solid
            this.shells.forEach(shell => {
                finalizedSolid.shells.push(shell);
            });
            // Perform final validation of the solid’s topology
            const validation = this.brepValidator.validateTopology(finalizedSolid);
            if (!validation) {
                throw new Error('Final topology validation failed for the solid');
            }
            return finalizedSolid;
        } catch (error) {
            throw new Error(`Finalization failed: ${ error.message }`);
        }
    }
}
class BrepKernel {
    constructor() {
        this.solids = [];
        this.shells = [];
        this.vertices = [];
        this.edges = [];
        this.faces = [];
        this.brepValidator = new BrepValidator();
        this.brepOperations = new BrepOperations();
        this.brepBuilder = new BrepBuilder();
    }
    create() {
        const newSolid = new Solid();
        this.solids.push(newSolid);
        return newSolid;
    }
    modify(entity) {
        // Input validation
        if (!(entity instanceof TopologyVertex) && !(entity instanceof TopologyEdge) && !(entity instanceof TopologyFace)) {
            throw new Error('Entity must be an instance of TopologyVertex, TopologyEdge, or TopologyFace');
        }
        try {
            // Validate the entity if it's part of a solid or shell
            const relatedEntity = this.findRelatedEntity(entity);
            if (!relatedEntity) {
                throw new Error('Entity is not part of any existing BREP structure');
            }
            // Modify the entity based on specific types
            if (entity instanceof TopologyVertex) {
                // Example modification logic for TopologyVertex
                entity.setPosition(entity.position.clone().add(new Vector3D(1, 0, 0)));
            } else // Example transformation
            if (entity instanceof TopologyEdge) {
                // Example modification logic for TopologyEdge
                entity.reverse();
            } else if (entity instanceof TopologyFace) {
                // Example modification logic for TopologyFace
                entity.addBound(new TopologyLoop());
            }
            // Example adding a new loop
            // Revalidate BREP structure after modification
            this.brepValidator.validateTopology(relatedEntity);
        } catch (error) {
            throw new Error(`Modification failed: ${ error.message }`);
        }
    }
    validate() {
        const validationResults = {
            isValid: true,
            errors: [],
            warnings: []
        };
        try {
            // Validate each solid in the kernel
            for (const solid of this.solids) {
                const shellValidationResult = this.brepValidator.validateTopology(solid);
                if (!shellValidationResult) {
                    validationResults.errors.push(`Solid validation failed for id: ${ solid.id }`);
                    validationResults.isValid = false;
                }
                solid.shells.forEach(shell => {
                    if (!this.brepValidator.validateTopology(shell)) {
                        validationResults.errors.push(`Shell validation failed for solid: ${ solid.id }`);
                        validationResults.isValid = false;
                    }
                });
            }
            // Check for empty structures
            if (this.solids.length === 0) {
                validationResults.errors.push('No solids have been created in the kernel');
                validationResults.isValid = false;
            }
            // Validate vertices and edges
            this.vertices.forEach(vertex => {
                if (!(vertex instanceof TopologyVertex)) {
                    validationResults.errors.push('Invalid vertex encountered in the kernel');
                    validationResults.isValid = false;
                }
            });
            this.edges.forEach(edge => {
                if (!(edge instanceof TopologyEdge)) {
                    validationResults.errors.push('Invalid edge encountered in the kernel');
                    validationResults.isValid = false;
                }
            });
            this.faces.forEach(face => {
                if (!(face instanceof TopologyFace)) {
                    validationResults.errors.push('Invalid face encountered in the kernel');
                    validationResults.isValid = false;
                }
            });
        } catch (error) {
            validationResults.errors.push(`Validation failed: ${ error.message }`);
            validationResults.isValid = false;
        }
        return Object.freeze(validationResults);
    }
    // ...existing methods
    import(file, format) {
        // Input validation
        if (typeof file !== 'string' || file.trim().length === 0) {
            throw new Error('File path must be a non-empty string');
        }
        if (format !== 'iges' && format !== 'step') {
            throw new Error('Invalid format specified. Supported formats are "iges" and "step"');
        }
        try {
            // Read file content (this part would need an actual file read, but here it's a placeholder)
            const data = this.readFileContent(file);
            // This method simulates file reading
            let geometries;
            if (format === 'iges') {
                const parsed = new FileIO().readIGES(data);
                geometries = parsed.geometries;
            } else if (format === 'step') {
                const parsed = new FileIO().readSTEP(data);
                geometries = parsed.geometries;
            }
            // Validate and store geometries in the kernel
            geometries.forEach(geometry => {
                if (geometry instanceof Solid) {
                    this.solids.push(geometry);
                } else if (geometry instanceof Shell) {
                    this.shells.push(geometry);
                } else if (geometry instanceof TopologyVertex) {
                    this.vertices.push(geometry);
                } else if (geometry instanceof TopologyEdge) {
                    this.edges.push(geometry);
                } else if (geometry instanceof TopologyFace) {
                    this.faces.push(geometry);
                } else {
                    throw new Error(`Unsupported geometry type: ${ geometry.constructor.name }`);
                }
            });
            // Validate the entire BREP structure
            const validationResults = this.validate();
            if (!validationResults.isValid) {
                throw new Error(`BREP import validation failed: ${ validationResults.errors.join(', ') }`);
            }
        } catch (error) {
            throw new Error(`Import failed: ${ error.message }`);
        }
    }
    // ...existing properties and methods
    export(format) {
        // Input validation
        if (format !== 'iges' && format !== 'step') {
            throw new Error('Invalid format specified. Supported formats are "iges" and "step"');
        }
        let geometries = [];
        try {
            // Gather all geometry entities from the kernel
            for (const solid of this.solids) {
                for (const shell of solid.shells) {
                    geometries.push(...shell.faces);
                    geometries.push(...shell.edges);
                    geometries.push(...shell.vertices);
                }
            }
            // Create an export handler
            const fileIO = new FileIO();
            let output;
            if (format === 'iges') {
                output = fileIO.writeIGES(geometries);
            } else if (format === 'step') {
                output = fileIO.writeSTEP(geometries);
            }
            // Return constructed file content
            return output;
        } catch (error) {
            throw new Error(`Export failed: ${ error.message }`);
        }
    }
    findRelatedEntity(entity) {
        // Helper method to find the related solid or shell containing the entity
        for (const solid of this.solids) {
            if (solid.shells.some(shell => shell.faces.includes(entity) || shell.edges.includes(entity))) {
                return solid;
            }
        }
        return null;
    }
    readFileContent(file) {
        // Placeholder for actual file reading logic
        return '';
    }
    optimize() {
        try {
            // Remove redundant vertices
            this.removeDuplicateVertices();
            // Simplify edges
            this.simplifyEdges();
            // Validate the final geometry state
            const validationResult = this.brepValidator.validateTopology(this);
            if (!validationResult) {
                throw new Error('BREP optimization validation failed');
            }
        } catch (error) {
            throw new Error(`Optimization failed: ${ error.message }`);
        }
    }
    // ... other methods ...
    removeDuplicateVertices() {
        try {
            // Map to store unique vertices based on position
            const uniqueVertices = new Map();
            const replacementMap = new Map();
            // First pass - identify unique vertices and build replacement map
            this.vertices.forEach(vertex => {
                const key = `${ vertex.position.x.toFixed(6) },${ vertex.position.y.toFixed(6) },${ vertex.position.z.toFixed(6) }`;
                if (!uniqueVertices.has(key)) {
                    uniqueVertices.set(key, vertex);
                } else {
                    // Map duplicate vertex to its replacement
                    replacementMap.set(vertex, uniqueVertices.get(key));
                }
            });
            // Second pass - update all topology references
            this.edges.forEach(edge => {
                if (replacementMap.has(edge.startVertex)) {
                    const replacement = replacementMap.get(edge.startVertex);
                    // Update edge reference
                    edge.startVertex = replacement;
                    // Update vertex's edge list
                    replacement.edges.push(edge);
                }
                if (replacementMap.has(edge.endVertex)) {
                    const replacement = replacementMap.get(edge.endVertex);
                    edge.endVertex = replacement;
                    replacement.edges.push(edge);
                }
            });
            // Update vertices collection with unique vertices only
            this.vertices = Array.from(uniqueVertices.values());
            // Clean up edge lists in remaining vertices
            this.vertices.forEach(vertex => {
                // Remove duplicates from edge list
                vertex.edges = Array.from(new Set(vertex.edges));
                // Remove any null or undefined references
                vertex.edges = vertex.edges.filter(edge => edge !== null && edge !== undefined);
            });
            // Validate the topology after removal
            const validationResult = this.brepValidator.validateTopology(this);
            if (!validationResult) {
                throw new Error('Topology validation failed after removing duplicate vertices');
            }
            return true;
        } catch (error) {
            throw new Error(`Failed to remove duplicate vertices: ${ error.message }`);
        }
    }
    simplifyEdges() {
        try {
            // Track edges that have been processed
            const processedEdges = new Set();
            const edgesToRemove = new Set();
            // Compare each edge with every other edge
            for (let i = 0; i < this.edges.length; i++) {
                const edge1 = this.edges[i];
                if (edgesToRemove.has(edge1))
                    continue;
                for (let j = i + 1; j < this.edges.length; j++) {
                    const edge2 = this.edges[j];
                    if (edgesToRemove.has(edge2))
                        continue;
                    // Check if edges are similar enough to merge
                    if (this.areEdgesSimilar(edge1, edge2)) {
                        // Merge the edges
                        const mergedEdge = this.brepOperations.mergeEdges(edge1, edge2);
                        if (mergedEdge) {
                            // Mark original edges for removal
                            edgesToRemove.add(edge1);
                            edgesToRemove.add(edge2);
                            // Add merged edge if it's not already present
                            if (!processedEdges.has(mergedEdge)) {
                                this.edges.push(mergedEdge);
                                processedEdges.add(mergedEdge);
                            }
                        }
                    }
                }
            }
            // Remove marked edges
            this.edges = this.edges.filter(edge => !edgesToRemove.has(edge));
            // Validate the topology after simplification
            const validationResult = this.brepValidator.validateTopology(this);
            if (!validationResult) {
                throw new Error('Topology validation failed after edge simplification');
            }
            return true;
        } catch (error) {
            throw new Error(`Edge simplification failed: ${ error.message }`);
        }
    }
    areEdgesSimilar(edge1, edge2) {
        try {
            // Check if edges have similar endpoints
            const startsSimilar = edge1.startVertex.position.subtract(edge2.startVertex.position).length() < defaultTolerance;
            const endsSimilar = edge1.endVertex.position.subtract(edge2.endVertex.position).length() < defaultTolerance;
            // Check for both direct and reversed similarity
            const directSimilar = startsSimilar && endsSimilar;
            const reversedSimilar = startsSimilar && edge1.startVertex.position.subtract(edge2.endVertex.position).length() < defaultTolerance;
            if (!directSimilar && !reversedSimilar)
                return false;
            // Check geometric similarity by sampling points along curves
            const numSamples = 10;
            const maxDeviation = defaultTolerance;
            for (let i = 0; i <= numSamples; i++) {
                const param = i / numSamples;
                const point1 = edge1.curve.evaluate(edge1.curve.domain.min + param * (edge1.curve.domain.max - edge1.curve.domain.min));
                const point2 = edge2.curve.evaluate(edge2.curve.domain.min + param * (edge2.curve.domain.max - edge2.curve.domain.min));
                if (point1.subtract(point2).length() > maxDeviation) {
                    return false;
                }
            }
            return true;
        } catch (error) {
            return false;
        }
    }
}
class BrepOperations {
    constructor() {
        this.brepValidator = new BrepValidator();
        this.brepBuilder = new BrepBuilder();
    }
    unite(solid1, solid2) {
        // Input validation
        if (!(solid1 instanceof Solid) || !(solid2 instanceof Solid)) {
            throw new Error('Both parameters must be instances of Solid');
        }
        try {
            // Create a new solid for the union result
            const unitedSolid = new Solid();
            // Combine shells from both solids
            solid1.shells.forEach(shell => {
                unitedSolid.shells.push(shell);
            });
            solid2.shells.forEach(shell => {
                unitedSolid.shells.push(shell);
            });
            // Validate the merged solid for topological integrity
            const validationResult = this.brepValidator.validateTopology(unitedSolid);
            if (!validationResult) {
                throw new Error('Topology validation failed after union operation');
            }
            return unitedSolid;
        } catch (error) {
            throw new Error(`Union operation failed: ${ error.message }`);
        }
    }
    subtract(solid1, solid2) {
        // Input validation
        if (!(solid1 instanceof Solid) || !(solid2 instanceof Solid)) {
            throw new Error('Both parameters must be instances of Solid');
        }
        try {
            // Create a new solid for the subtraction result
            const resultSolid = new Solid();
            // Iterate over faces of the first solid
            for (const face1 of solid1.shells.flatMap(shell => shell.faces)) {
                // Iterate over faces of the second solid
                for (const face2 of solid2.shells.flatMap(shell => shell.faces)) {
                    // Check for intersection between the two faces
                    const intersectionResult = this.calculateFaceIntersection(face1, face2);
                    // If there is an intersection, handle the subtraction
                    if (intersectionResult.points.length > 0) {
                        this.handleFaceSubtraction(resultSolid, face1, intersectionResult);
                    } else {
                        // If no intersection, add face1 to result solid
                        resultSolid.addFace(face1);
                    }
                }
            }
            // Validate the resulting solid's topology
            const validationResult = this.brepValidator.validateTopology(resultSolid);
            if (!validationResult) {
                throw new Error('Topology validation failed after subtraction operation');
            }
            return resultSolid;
        } catch (error) {
            throw new Error(`Subtraction operation failed: ${ error.message }`);
        }
    }
    intersect(solid1, solid2) {
        // Input validation
        if (!(solid1 instanceof Solid) || !(solid2 instanceof Solid)) {
            throw new Error('Both parameters must be instances of Solid');
        }
        try {
            // Create a new solid for the intersection result
            const intersectionSolid = new Solid();
            // Iterate over faces of the first solid
            for (const face1 of solid1.shells.flatMap(shell => shell.faces)) {
                // Iterate over faces of the second solid
                for (const face2 of solid2.shells.flatMap(shell => shell.faces)) {
                    // Check for intersection between the two faces
                    const intersectionResult = this.calculateFaceIntersection(face1, face2);
                    // If there is an intersection, handle the creation of new faces
                    if (intersectionResult.points.length > 0) {
                        this.handleFaceIntersection(intersectionSolid, face1, intersectionResult);
                    }
                }
            }
            // Validate the resulting solid's topology
            const validationResult = this.brepValidator.validateTopology(intersectionSolid);
            if (!validationResult) {
                throw new Error('Topology validation failed after intersection operation');
            }
            return intersectionSolid;
        } catch (error) {
            throw new Error(`Intersection operation failed: ${ error.message }`);
        }
    }
    imprint(face) {
        // Input validation
        if (!(face instanceof TopologyFace)) {
            throw new Error('Parameter must be a TopologyFace instance');
        }
        try {
            // Create a new face that will be the result of the imprint operation
            const imprintedFace = new TopologyFace();
            const newEdges = [];
            // Iterate over the trim curves of the face
            for (const curve of face.edges) {
                // Calculate intersections with the face
                const intersections = curve.getIntersections(face);
                if (intersections.length > 0) {
                    for (const intersection of intersections) {
                        // Create new edge segments based on intersections
                        const newEdge = new TopologyEdge();
                        newEdge.curve = intersection.curve;
                        // Store the intersection curve
                        newEdges.push(newEdge);
                    }
                }
            }
            // Add the new edges to the imprinted face
            imprintedFace.edges = newEdges;
            // Validate the resulting face
            imprintedFace.validateBounds();
            // Return the new imprinted face
            return imprintedFace;
        } catch (error) {
            throw new Error(`Imprint operation failed: ${ error.message }`);
        }
    }
    split() {
        // Input validation
        if (this.shells.length === 0) {
            throw new Error('No shells available for splitting');
        }
        try {
            // Create an array to store the new solids created from the split
            const newSolids = [];
            // Iterate through each shell in the current solid
            for (const shell of this.shells) {
                // Validate shell for splitting
                shell.validate();
                // Create a new shell for split results
                const newShell = new Shell();
                // Split each face in the shell based on specified criteria
                shell.faces.forEach(face => {
                    const splitFaces = face.split(edge => {
                        // Criteria for splitting, e.g., based on edge lengths or adjacency
                        return edge.length() > defaultTolerance;
                    });
                    // Placeholder criteria
                    // Add the newly created faces to the new shell
                    splitFaces.forEach(newFace => newShell.addFace(newFace));
                });
                // Add the new shell to the array of new solids
                newSolids.push(newShell);
            }
            // Return the new solids created from the split operation
            return newSolids;
        } catch (error) {
            throw new Error(`Split operation failed: ${ error.message }`);
        }
    }
    merge(solid1, solid2) {
        // Input validation
        if (!(solid1 instanceof Solid) || !(solid2 instanceof Solid)) {
            throw new Error('Both parameters must be instances of Solid');
        }
        try {
            // Create a new solid for the merged result
            const mergedSolid = new Solid();
            // Combine shells from both solids
            for (const shell1 of solid1.shells) {
                mergedSolid.shells.push(shell1.clone());
            }
            for (const shell2 of solid2.shells) {
                mergedSolid.shells.push(shell2.clone());
            }
            // Check for and merge overlapping edges and faces
            const facePairs = [];
            for (const face1 of mergedSolid.shells.flatMap(shell => shell.faces)) {
                for (const face2 of mergedSolid.shells.flatMap(shell => shell.faces)) {
                    if (face1 !== face2) {
                        const intersectionResult = this.calculateFaceIntersection(face1, face2);
                        if (intersectionResult.points.length > 0) {
                            facePairs.push({
                                face1,
                                face2,
                                intersectionResult
                            });
                        }
                    }
                }
            }
            for (const {face1, face2, intersectionResult} of facePairs) {
                this.handleFaceSubtraction(mergedSolid, face1, intersectionResult);
                this.handleFaceIntersection(mergedSolid, face1, intersectionResult);
            }
            // Validate the final merged solid's topology
            const validationResult = this.brepValidator.validateTopology(mergedSolid);
            if (!validationResult) {
                throw new Error('Topology validation failed after merge operation');
            }
            return mergedSolid;
        } catch (error) {
            throw new Error(`Merge operation failed: ${ error.message }`);
        }
    }
    calculateFaceIntersection(face1, face2) {
        // Placeholder for actual intersection calculation logic
        return { points: [] };
    }
    handleFaceSubtraction(resultSolid, face1, intersectionResult) {
        try {
            // Validate inputs
            if (!(resultSolid instanceof Solid)) {
                throw new Error('Result must be a Solid instance');
            }
            if (!(face1 instanceof TopologyFace)) {
                throw new Error('Face1 must be a TopologyFace instance');
            }
            if (!intersectionResult || !intersectionResult.points || intersectionResult.points.length === 0) {
                throw new Error('Intersection result must contain points');
            }
            // Create new vertices from intersection points
            const intersectionVertices = intersectionResult.points.map(point => new TopologyVertex({ position: point }));
            // Create edges along intersection curve
            const intersectionEdges = [];
            for (let i = 0; i < intersectionVertices.length - 1; i++) {
                const edge = new TopologyEdge();
                const controlPoints = [
                    new ControlPoint(intersectionVertices[i].position.x, intersectionVertices[i].position.y, intersectionVertices[i].position.z),
                    new ControlPoint(intersectionVertices[i + 1].position.x, intersectionVertices[i + 1].position.y, intersectionVertices[i + 1].position.z)
                ];
                // Create linear NURBS curve for intersection edge
                const knotVector = new KnotVector([
                    0,
                    0,
                    1,
                    1
                ]);
                edge.curve = new NurbsCurve(controlPoints, knotVector, 1);
                edge.startVertex = intersectionVertices[i];
                edge.endVertex = intersectionVertices[i + 1];
                intersectionEdges.push(edge);
            }
            // Split face using intersection edges
            const newFaces = [];
            let currentFace = face1.clone();
            for (const edge of intersectionEdges) {
                // Find parameter for splitting
                const closestPoint = currentFace.surface.closestPoint(edge.startVertex.position);
                if (closestPoint) {
                    // Split current face
                    const splitResult = currentFace.split(closestPoint.parameters);
                    // Keep the part that's not intersecting
                    for (const splitFace of splitResult) {
                        // Check if this part is outside the intersection
                        const centroid = splitFace.calculateCentroid();
                        if (!this.isPointInIntersection(centroid, intersectionResult)) {
                            newFaces.push(splitFace);
                        }
                    }
                    // Update current face for next iteration
                    currentFace = splitResult[splitResult.length - 1];
                }
            }
            // Add non-intersecting faces to result solid
            for (const newFace of newFaces) {
                resultSolid.addFace(newFace);
            }
            // Validate the modified solid
            if (!this.brepValidator.validateTopology(resultSolid)) {
                throw new Error('Invalid topology after face subtraction');
            }
            return true;
        } catch (error) {
            throw new Error(`Failed to handle face subtraction: ${ error.message }`);
        }
    }
    handleFaceIntersection(resultSolid, face1, intersectionResult) {
        // Input validation
        if (!(resultSolid instanceof Solid)) {
            throw new Error('Result must be a Solid instance');
        }
        if (!(face1 instanceof TopologyFace)) {
            throw new Error('Face1 must be a TopologyFace instance');
        }
        if (!intersectionResult || !intersectionResult.points || intersectionResult.points.length === 0) {
            throw new Error('Intersection result must contain points');
        }
        try {
            // Create a new face based on the intersection curve
            const intersectionFace = new TopologyFace();
            // Create vertices for intersection points
            const vertices = intersectionResult.points.map(point => new TopologyVertex({ position: point }));
            // Create edges connecting the intersection points
            for (let i = 0; i < vertices.length - 1; i++) {
                const edge = new TopologyEdge();
                // Create a NURBS curve between consecutive points
                const controlPoints = [
                    new ControlPoint(vertices[i].position.x, vertices[i].position.y, vertices[i].position.z),
                    new ControlPoint(vertices[i + 1].position.x, vertices[i + 1].position.y, vertices[i + 1].position.z)
                ];
                const knotVector = new KnotVector([
                    0,
                    0,
                    1,
                    1
                ]);
                // Degree 1 curve
                edge.curve = new NurbsCurve(controlPoints, knotVector, 1);
                // Set edge vertices
                edge.startVertex = vertices[i];
                edge.endVertex = vertices[i + 1];
                // Update vertex references
                vertices[i].edges.push(edge);
                vertices[i + 1].edges.push(edge);
                // Add edge to face
                intersectionFace.edges.push(edge);
            }
            // Create closing edge if needed (for closed intersection curves)
            if (vertices.length > 2 && vertices[0].position.subtract(vertices[vertices.length - 1].position).length() < defaultTolerance) {
                const closingEdge = new TopologyEdge();
                const controlPoints = [
                    new ControlPoint(vertices[vertices.length - 1].position.x, vertices[vertices.length - 1].position.y, vertices[vertices.length - 1].position.z),
                    new ControlPoint(vertices[0].position.x, vertices[0].position.y, vertices[0].position.z)
                ];
                const knotVector = new KnotVector([
                    0,
                    0,
                    1,
                    1
                ]);
                closingEdge.curve = new NurbsCurve(controlPoints, knotVector, 1);
                closingEdge.startVertex = vertices[vertices.length - 1];
                closingEdge.endVertex = vertices[0];
                vertices[vertices.length - 1].edges.push(closingEdge);
                vertices[0].edges.push(closingEdge);
                intersectionFace.edges.push(closingEdge);
            }
            // Validate face edges
            intersectionFace.validateBounds();
            // Add the intersection face to the result solid
            resultSolid.addFace(intersectionFace);
        } catch (error) {
            throw new Error(`Failed to handle face intersection: ${ error.message }`);
        }
    }
    mergeEdges(edge1, edge2) {
        // Input validation
        if (!(edge1 instanceof TopologyEdge) || !(edge2 instanceof TopologyEdge)) {
            throw new Error('Both parameters must be TopologyEdge instances');
        }
        try {
            // Check if edges can be merged (share common vertex)
            const commonVertex = this.findCommonVertex(edge1, edge2);
            if (!commonVertex) {
                throw new Error('Edges must share a common vertex to merge');
            }
            // Create new merged edge
            const mergedEdge = new TopologyEdge();
            // Determine merge direction and set vertices
            if (edge1.endVertex === commonVertex && edge2.startVertex === commonVertex) {
                mergedEdge.startVertex = edge1.startVertex;
                mergedEdge.endVertex = edge2.endVertex;
                // Join curves in correct order
                mergedEdge.curve = edge1.curve.join(edge2.curve);
            } else if (edge1.startVertex === commonVertex && edge2.endVertex === commonVertex) {
                mergedEdge.startVertex = edge2.startVertex;
                mergedEdge.endVertex = edge1.endVertex;
                mergedEdge.curve = edge2.curve.join(edge1.curve);
            } else if (edge1.startVertex === commonVertex && edge2.startVertex === commonVertex) {
                // Need to reverse one edge
                edge2.reverse();
                mergedEdge.startVertex = edge2.endVertex;
                mergedEdge.endVertex = edge1.endVertex;
                mergedEdge.curve = edge2.curve.join(edge1.curve);
            } else if (edge1.endVertex === commonVertex && edge2.endVertex === commonVertex) {
                // Need to reverse one edge
                edge2.reverse();
                mergedEdge.startVertex = edge1.startVertex;
                mergedEdge.endVertex = edge2.startVertex;
                mergedEdge.curve = edge1.curve.join(edge2.curve);
            }
            // Update vertex references
            mergedEdge.startVertex.edges = mergedEdge.startVertex.edges.filter(e => e !== edge1 && e !== edge2);
            mergedEdge.startVertex.edges.push(mergedEdge);
            mergedEdge.endVertex.edges = mergedEdge.endVertex.edges.filter(e => e !== edge1 && e !== edge2);
            mergedEdge.endVertex.edges.push(mergedEdge);
            // Merge face references
            mergedEdge.faces = Array.from(new Set([
                ...edge1.faces,
                ...edge2.faces
            ]));
            // Update face references
            mergedEdge.faces.forEach(face => {
                const index1 = face.edges.indexOf(edge1);
                const index2 = face.edges.indexOf(edge2);
                if (index1 !== -1)
                    face.edges[index1] = mergedEdge;
                if (index2 !== -1)
                    face.edges.splice(index2, 1);
            });
            // Validate merged edge
            this.validateMergedEdge(mergedEdge);
            // Invalidate original edges
            edge1.curve = null;
            edge1.startVertex = null;
            edge1.endVertex = null;
            edge1.faces = [];
            edge2.curve = null;
            edge2.startVertex = null;
            edge2.endVertex = null;
            edge2.faces = [];
            return mergedEdge;
        } catch (error) {
            throw new Error(`Edge merge operation failed: ${ error.message }`);
        }
    }
    isPointInIntersection(point, intersectionResult) {
        // Check if a point lies within the intersection region
        for (const intersectionPoint of intersectionResult.points) {
            if (point.subtract(intersectionPoint).length() < defaultTolerance) {
                return true;
            }
        }
        return false;
    }
    findCommonVertex(edge1, edge2) {
        // Input validation
        if (!(edge1 instanceof TopologyEdge) || !(edge2 instanceof TopologyEdge)) {
            throw new Error('Both parameters must be TopologyEdge instances');
        }
        try {
            // Compare all possible vertex combinations
            if (edge1.startVertex === edge2.startVertex)
                return edge1.startVertex;
            if (edge1.startVertex === edge2.endVertex)
                return edge1.startVertex;
            if (edge1.endVertex === edge2.startVertex)
                return edge1.endVertex;
            if (edge1.endVertex === edge2.endVertex)
                return edge1.endVertex;
            // If no exact match, check for vertices within tolerance
            const tolerance = Math.min(edge1.startVertex.tolerance, edge1.endVertex.tolerance, edge2.startVertex.tolerance, edge2.endVertex.tolerance);
            // Check all combinations with tolerance
            const vertices = [
                [
                    edge1.startVertex,
                    edge2.startVertex
                ],
                [
                    edge1.startVertex,
                    edge2.endVertex
                ],
                [
                    edge1.endVertex,
                    edge2.startVertex
                ],
                [
                    edge1.endVertex,
                    edge2.endVertex
                ]
            ];
            for (const [v1, v2] of vertices) {
                if (v1.position.subtract(v2.position).length() <= tolerance) {
                    // If vertices are within tolerance, merge them and return the result
                    return v1.merge(v2);
                }
            }
            // No common vertex found
            return null;
        } catch (error) {
            throw new Error(`Failed to find common vertex: ${ error.message }`);
        }
    }
}
class BrepValidator {
    constructor() {
        this.errors = [];
        this.warnings = [];
    }
    validateTopology(shell) {
        // Input validation
        if (!(shell instanceof Shell)) {
            throw new Error('Parameter must be a Shell instance');
        }
        // Check if shell has at least one face
        if (shell.faces.length === 0) {
            this.errors.push('Shell must contain at least one face');
            return false;
        }
        // Validate each face in the shell
        for (const face of shell.faces) {
            // Check for edge belonging to the face
            if (face.edges.length === 0) {
                this.errors.push('Face must have at least one edge');
                return false;
            }
            // Validate edges consistency
            for (const edge of face.edges) {
                if (!edge.faces.includes(face)) {
                    this.errors.push('Edge must belong to the face');
                    return false;
                }
            }
        }
        // Validate edge connectivity across faces
        const edgeMap = new Map();
        for (const face of shell.faces) {
            for (const edge of face.edges) {
                const count = edgeMap.get(edge) || 0;
                edgeMap.set(edge, count + 1);
            }
        }
        for (const [edge, count] of edgeMap) {
            if (count !== 2) {
                this.errors.push('Non-manifold edge detected: edge must belong to exactly two faces');
                return false;
            }
        }
        // If no errors were found, return true
        return this.errors.length === 0;
    }
    validateGeometry() {
        // Check if all vertices have valid positions connected by edges
        const vertexSet = new Set(this.vertices);
        for (const vertex of vertexSet) {
            if (!(vertex instanceof TopologyVertex)) {
                this.errors.push('Invalid vertex in geometry');
                continue;
            }
            if (!vertex.position) {
                this.errors.push(`Vertex ${ vertex.id } has no position defined`);
                continue;
            }
            const edgeCount = vertex.edges.length;
            if (edgeCount === 0) {
                this.errors.push(`Vertex ${ vertex.id } is not connected to any edges`);
            }
        }
        // Check if all edges are connected to valid vertices
        const edgeSet = new Set(this.edges);
        for (const edge of edgeSet) {
            if (!(edge instanceof TopologyEdge)) {
                this.errors.push('Invalid edge in geometry');
                continue;
            }
            if (!edge.startVertex || !edge.endVertex) {
                this.errors.push(`Edge is not properly connected - missing vertices`);
                continue;
            }
            if (!vertexSet.has(edge.startVertex) || !vertexSet.has(edge.endVertex)) {
                this.errors.push(`Edge ${ edge.startVertex.id }-${ edge.endVertex.id } connects invalid vertices`);
            }
        }
        // Check if all faces are valid
        const faceSet = new Set(this.faces);
        for (const face of faceSet) {
            if (!(face instanceof TopologyFace)) {
                this.errors.push('Invalid face in geometry');
                continue;
            }
            if (face.edges.length < 3) {
                this.errors.push(`Face with edges ${ face.edges.length } is not valid - must have at least three edges`);
            }
            const edgeCount = face.edges.length;
            for (const edge of face.edges) {
                if (!edgeSet.has(edge)) {
                    this.errors.push(`Face ${ face } references an invalid edge`);
                }
            }
        }
        return this.errors.length === 0;
    }
    validateOrientation(shell) {
        // Input validation
        if (!(shell instanceof Shell)) {
            throw new Error('Parameter must be a Shell instance');
        }
        // Check if shell has at least one face
        if (shell.faces.length === 0) {
            this.errors.push('Shell must contain at least one face to validate orientation');
            return false;
        }
        const faceNormals = new Map();
        // Validate orientation for each face
        for (const face of shell.faces) {
            // Ensure the face has associated edges
            if (face.edges.length === 0) {
                this.errors.push('Face has no edges, cannot validate orientation');
                return false;
            }
            // Calculate normal for the face using cross product of edges
            const edge1 = face.edges[0];
            const edge2 = face.edges[1];
            const startPoint = edge1.startVertex.position;
            const endPoint1 = edge1.endVertex.position;
            const endPoint2 = edge2.endVertex.position;
            const normal = endPoint1.subtract(startPoint).cross(endPoint2.subtract(startPoint)).normalize();
            // Check if normal is already recorded
            const normalKey = `${ normal.x.toFixed(6) },${ normal.y.toFixed(6) },${ normal.z.toFixed(6) }`;
            if (!faceNormals.has(normalKey)) {
                faceNormals.set(normalKey, face);
            } else {
                // Check if the current face has the same orientation as previously recorded face
                const existingFace = faceNormals.get(normalKey);
                if (!this.areNormalsConsistent(normal, existingFace)) {
                    this.errors.push('Inconsistent face orientations detected in the shell');
                    return false;
                }
            }
        }
        // All checks passed
        return true;
    }
    healTopology(shell) {
        // Input validation
        if (!(shell instanceof Shell)) {
            throw new Error('Parameter must be a Shell instance');
        }
        // Attempt to heal edges within tolerance
        for (const face of shell.faces) {
            for (let i = 0; i < face.edges.length; i++) {
                const edge = face.edges[i];
                const nextEdge = face.edges[(i + 1) % face.edges.length];
                // Check for gaps between this edge and the next
                const endVertex = edge.endVertex.position;
                const startVertex = nextEdge.startVertex.position;
                if (endVertex.distanceTo(startVertex) > defaultTolerance) {
                    // Calculate midpoint
                    const midPoint = new Vector3D((endVertex.x + startVertex.x) / 2, (endVertex.y + startVertex.y) / 2, (endVertex.z + startVertex.z) / 2);
                    // Attempt to set position of start vertex of next edge
                    try {
                        nextEdge.startVertex.setPosition(midPoint);
                    } catch (error) {
                        this.warnings.push(`Failed to heal edge between vertices ${ edge.startVertex.id } and ${ nextEdge.endVertex.id }: ${ error.message }`);
                    }
                }
            }
        }
    }
    areNormalsConsistent(normal, existingFace) {
        // Check for consistent normal direction
        return Math.abs(normal.dot(existingFace.getNormal())) > 0.99999;
    }
}
class ControlPoint {
    constructor(x = 0, y = 0, z = 0, w = 1) {
        this.validateCoordinates(x, y, z, w);
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
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
    // Stub for any missing methods
    fromHomogeneous(coords) {
        // Input validation
        if (!Array.isArray(coords) || coords.length !== 4) {
            throw new Error('Input must be an array of four elements [x, y, z, w]');
        }
        const [x, y, z, w] = coords.map(coord => {
            if (typeof coord !== 'number' || !Number.isFinite(coord)) {
                throw new Error('All coordinates must be finite numbers');
            }
            return coord;
        });
        // Validate weight (w) is non-zero
        if (Math.abs(w) < Number.EPSILON) {
            throw new Error('Weight (w) must be non-zero for valid homogeneous coordinates');
        }
        // Assign the converted coordinates
        return new ControlPoint(x * w, y * w, z * w, w);
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
    // ...existing methods...
    clone() {
        try {
            return new ControlPoint(this.x === 0 ? 0 : this.x, this.y === 0 ? 0 : this.y, this.z === 0 ? 0 : this.z, this.w === 0 ? 1 : this.w);
        } catch (error) {
            throw new Error(`Failed to clone control point: ${ error.message }`);
        }
    }
    validateCoordinates(x, y, z, w) {
        if (![
                x,
                y,
                z,
                w
            ].every(Number.isFinite)) {
            throw new Error('All coordinates and weight must be finite numbers');
        }
        if (Math.abs(w) < Number.EPSILON) {
            throw new Error('Weight (w) must be non-zero for valid homogeneous coordinates');
        }
    }
}
class FileIO {
    constructor(options = {}) {
        // Default options
        const defaultOptions = {
            encoding: 'utf-8',
            precision: 6,
            maxFileSize: 1024 * 1024 * 50,
            // 50MB
            supportedFormats: {
                iges: [
                    '.igs',
                    '.iges'
                ],
                step: [
                    '.stp',
                    '.step',
                    '.p21'
                ]
            },
            validateInput: true,
            preserveMetadata: true
        };
        // Merge provided options with defaults
        this.options = {
            ...defaultOptions,
            ...options
        };
        // Initialize internal state
        this.lastError = null;
        this.currentFile = null;
        this.geometryCache = new Map();
        // File format handlers
        this.formatHandlers = {
            iges: {
                read: null,
                // Will be initialized when needed
                write: null
            },
            // Will be initialized when needed
            step: {
                read: null,
                // Will be initialized when needed
                write: null
            }
        };
        // Will be initialized when needed
        // Validate options
        this.validateOptions();
        // Freeze critical properties to prevent modification
        Object.freeze(this.options.supportedFormats);
    }
    readIGES(data) {
        try {
            // Input validation
            if (!data || typeof data !== 'string') {
                throw new Error('Input data must be a non-empty string');
            }
            // Initialize storage for parsed entities
            const entities = [];
            const parameters = new Map();
            const geometries = [];
            // Split data into sections based on IGES format
            const lines = data.split('\n');
            const sections = {
                start: [],
                global: [],
                directory: [],
                parameter: [],
                terminate: null
            };
            // Parse sections based on last character of each line
            for (const line of lines) {
                if (line.length < 73)
                    continue;
                const marker = line.charAt(72);
                switch (marker) {
                case 'S':
                    sections.start.push(line);
                    break;
                case 'G':
                    sections.global.push(line);
                    break;
                case 'D':
                    sections.directory.push(line);
                    break;
                case 'P':
                    sections.parameter.push(line);
                    break;
                case 'T':
                    sections.terminate = line;
                    break;
                }
            }
            // Parse global section for units and other settings
            const globalData = {
                modelSpaceScale: 1,
                unitFlag: 1,
                // Default to millimeters
                maxCoordValue: 0,
                author: '',
                organization: ''
            };
            if (sections.global.length > 0) {
                const globalSection = sections.global.join('').substring(0, 72);
                globalData.modelSpaceScale = parseFloat(globalSection.substring(48, 56)) || 1;
                globalData.unitFlag = parseInt(globalSection.substring(12, 24)) || 1;
                globalData.maxCoordValue = parseFloat(globalSection.substring(24, 36)) || 0;
            }
            // Parse directory entries
            for (let i = 0; i < sections.directory.length; i += 2) {
                const entry1 = sections.directory[i];
                const entry2 = i + 1 < sections.directory.length ? sections.directory[i + 1] : '';
                if (entry1.length >= 72 && entry2.length >= 72) {
                    const entityType = parseInt(entry1.substring(0, 8));
                    const parameterStart = parseInt(entry1.substring(8, 16));
                    const structure = parseInt(entry1.substring(16, 24));
                    const lineFont = parseInt(entry1.substring(24, 32));
                    const level = parseInt(entry1.substring(32, 40));
                    const view = parseInt(entry1.substring(40, 48));
                    const matrix = parseInt(entry1.substring(48, 56));
                    const label = parseInt(entry1.substring(56, 64));
                    const status = entry1.substring(64, 72);
                    const color = parseInt(entry2.substring(8, 16));
                    const parameterLineCount = parseInt(entry2.substring(16, 24));
                    const form = parseInt(entry2.substring(24, 32));
                    entities.push({
                        entityType,
                        parameterStart,
                        structure,
                        lineFont,
                        level,
                        view,
                        matrix,
                        label,
                        status,
                        color,
                        parameterLineCount,
                        form
                    });
                }
            }
            // Parse parameter data
            let currentParam = '';
            for (const line of sections.parameter) {
                // Remove section marker and sequence number
                const paramData = line.substring(0, 64);
                currentParam += paramData;
                if (line.charAt(71) === '1') {
                    // End of parameter marker
                    const paramValues = currentParam.split(',').map(p => p.trim()).filter(p => p.length > 0);
                    parameters.set(parameters.size + 1, paramValues);
                    currentParam = '';
                }
            }
            // Convert entities to geometries
            for (const entity of entities) {
                const paramData = parameters.get(entity.parameterStart);
                if (!paramData)
                    continue;
                switch (entity.entityType) {
                case 126:
                    // NURBS Curve
                    geometries.push(this.parseNurbsCurve(paramData, globalData));
                    break;
                case 128:
                    // NURBS Surface
                    geometries.push(this.parseNurbsSurface(paramData, globalData));
                    break;
                case 144:
                    // Trimmed Surface
                    geometries.push(this.parseTrimmedSurface(paramData, globalData, entities, parameters));
                    break;
                }
            }
            // Return parsed geometries
            return Object.freeze({
                geometries: Object.freeze(geometries),
                globalData: Object.freeze(globalData),
                success: true,
                metadata: Object.freeze({
                    fileType: 'IGES',
                    entityCount: entities.length,
                    geometryCount: geometries.length,
                    processedDate: new Date()
                })
            });
        } catch (error) {
            throw new Error(`IGES parsing failed: ${ error.message }`);
        }
    }
    writeIGES(geometries) {
        try {
            // Input validation
            if (!Array.isArray(geometries)) {
                throw new Error('Input must be an array of geometries');
            }
            const startSection = [];
            const globalSection = [];
            const directorySection = [];
            const parameterSection = [];
            let parameterLineCount = 0;
            let entityNumber = 1;
            // Add start section header
            startSection.push(this.formatIGESLine('                                                                        ', 'S', 1));
            // Add global section
            const globalParams = [
                '1H,',
                // Parameter delimiter
                '1H;',
                // Record delimiter
                '4HIGES',
                // File name
                '4H1.0,',
                // System ID
                '32,',
                // Max line width
                '38,',
                // Max line width for parameter data
                '6,',
                // Sender's system ID
                '308,',
                // Product ID from sender
                '15,',
                // File version
                '4,',
                // Draft standard version
                ',',
                // Creation date/time
                ',',
                // Min resolution
                '1.0,',
                // Max coordinate value
                '1.0,',
                // Author
                ',',
                // Organization
                '11,',
                // Version flag
                '0;'
            ].join('');
            globalSection.push(this.formatIGESLine(globalParams, 'G', 1));
            // Process each geometry
            for (const geometry of geometries) {
                let entityData;
                if (geometry instanceof NurbsCurve) {
                    entityData = this.writeNurbsCurve(geometry, entityNumber, parameterLineCount + 1);
                } else if (geometry instanceof NurbsSurface) {
                    entityData = this.writeNurbsSurface(geometry, entityNumber, parameterLineCount + 1);
                } else if (geometry instanceof TrimmedSurface) {
                    // Handle trimmed surface (would need additional implementation)
                    continue;
                } else {
                    throw new Error(`Unsupported geometry type: ${ geometry.constructor.name }`);
                }
                // Add entity data to respective sections
                directorySection.push(...entityData.directory);
                parameterSection.push(...entityData.parameter);
                // Update counters
                parameterLineCount += entityData.parameter.length;
                entityNumber++;
            }
            // Create terminate section
            const terminateSection = this.formatIGESLine(`S${ String(startSection.length).padStart(7) }G${ String(globalSection.length).padStart(7) }` + `D${ String(directorySection.length).padStart(7) }P${ String(parameterSection.length).padStart(7) }T${ String(1).padStart(7) }`, 'T', 1);
            // Combine all sections
            const igesFile = [
                ...startSection,
                ...globalSection,
                ...directorySection,
                ...parameterSection,
                terminateSection
            ].join('\n');
            // Return IGES file content
            return igesFile;
        } catch (error) {
            throw new Error(`IGES file generation failed: ${ error.message }`);
        }
    }
    readSTEP(data) {
        try {
            // Input validation
            if (!data || typeof data !== 'string') {
                throw new Error('Input data must be a non-empty string');
            }
            // Initialize storage for parsed entities
            const entities = new Map();
            const geometries = [];
            // Split data into lines and remove comments
            const lines = data.split('\n').map(line => line.split('/*')[0].trim()).filter(line => line.length > 0);
            // Parse header section
            const headerData = {
                fileDescription: '',
                fileName: '',
                fileSchema: '',
                timeStamp: null
            };
            let currentSection = '';
            let currentEntity = null;
            // Process each line
            for (let i = 0; i < lines.length; i++) {
                const line = lines[i];
                // Check for section markers
                if (line.startsWith('HEADER;')) {
                    currentSection = 'HEADER';
                    continue;
                } else if (line.startsWith('DATA;')) {
                    currentSection = 'DATA';
                    continue;
                } else if (line.startsWith('ENDSEC;')) {
                    currentSection = '';
                    continue;
                }
                // Process header section
                if (currentSection === 'HEADER') {
                    if (line.includes('FILE_DESCRIPTION')) {
                        headerData.fileDescription = this.parseHeaderAttribute(line);
                    } else if (line.includes('FILE_NAME')) {
                        headerData.fileName = this.parseHeaderAttribute(line);
                        headerData.timeStamp = this.parseTimeStamp(line);
                    } else if (line.includes('FILE_SCHEMA')) {
                        headerData.fileSchema = this.parseHeaderAttribute(line);
                    }
                }
                // Process data section
                if (currentSection === 'DATA') {
                    if (line.startsWith('#')) {
                        // New entity definition
                        if (currentEntity) {
                            this.processEntity(currentEntity, entities);
                        }
                        currentEntity = {
                            id: this.parseEntityId(line),
                            type: this.parseEntityType(line),
                            data: line
                        };
                    } else if (currentEntity && !line.endsWith(';')) {
                        // Continuation of current entity
                        currentEntity.data += line;
                    } else if (currentEntity) {
                        // End of current entity
                        currentEntity.data += line;
                        this.processEntity(currentEntity, entities);
                        currentEntity = null;
                    }
                }
            }
            // Process any remaining entity
            if (currentEntity) {
                this.processEntity(currentEntity, entities);
            }
            // Convert entities to geometries
            for (const [id, entity] of entities) {
                const geometry = this.convertToGeometry(entity, entities);
                if (geometry) {
                    geometries.push(geometry);
                }
            }
            // Return parsed data
            return Object.freeze({
                geometries: Object.freeze(geometries),
                header: Object.freeze(headerData),
                success: true,
                metadata: Object.freeze({
                    fileType: 'STEP',
                    entityCount: entities.size,
                    geometryCount: geometries.length,
                    processedDate: new Date()
                })
            });
        } catch (error) {
            throw new Error(`STEP parsing failed: ${ error.message }`);
        }
    }
    writeSTEP(geometries, options = {}) {
        // Input validation
        if (!Array.isArray(geometries)) {
            throw new Error('Input must be an array of geometries');
        }
        // Default options
        const defaultOptions = {
            schema: 'AP214',
            precision: 6,
            includeMetadata: true
        };
        const settings = {
            ...defaultOptions,
            ...options
        };
        try {
            // Initialize sections
            const headerSection = [];
            const dataSection = [];
            let entityCounter = 1;
            // Generate header section
            headerSection.push('ISO-10303-21;');
            headerSection.push('HEADER;');
            headerSection.push('FILE_DESCRIPTION((\'' + settings.schema + ' CONFIGURATION CONTROLLED 3D DESIGNS\'), \'2;1\');');
            // Add file name and timestamp
            const timestamp = new Date().toISOString();
            const fileName = `FILE_NAME('${ options.fileName || 'export.stp' }','${ timestamp }',('AUTHOR'),('ORGANIZATION'),'','','');`;
            headerSection.push(fileName);
            // Add schema declaration
            headerSection.push(`FILE_SCHEMA(('${ settings.schema }'));`);
            headerSection.push('ENDSEC;');
            // Start data section
            dataSection.push('DATA;');
            // Process each geometry
            for (const geometry of geometries) {
                if (geometry instanceof NurbsCurve) {
                    const curveEntities = this.writeNurbsCurveSTEP(geometry, entityCounter, settings);
                    dataSection.push(...curveEntities.data);
                    entityCounter = curveEntities.nextId;
                } else if (geometry instanceof NurbsSurface) {
                    const surfaceEntities = this.writeNurbsSurfaceSTEP(geometry, entityCounter, settings);
                    dataSection.push(...surfaceEntities.data);
                    entityCounter = surfaceEntities.nextId;
                } else if (geometry instanceof TrimmedSurface) {
                    const trimmedEntities = this.writeTrimmedSurfaceSTEP(geometry, entityCounter, settings);
                    dataSection.push(...trimmedEntities.data);
                    entityCounter = trimmedEntities.nextId;
                }
            }
            // End sections
            dataSection.push('ENDSEC;');
            dataSection.push('END-ISO-10303-21;');
            // Combine all sections
            const fileContent = [
                ...headerSection,
                ...dataSection
            ].join('\n');
            return fileContent;
        } catch (error) {
            throw new Error(`STEP file generation failed: ${ error.message }`);
        }
    }
    validateOptions() {
        // Validate encoding
        if (typeof this.options.encoding !== 'string' || this.options.encoding.trim().length === 0) {
            throw new Error('Invalid encoding specified');
        }
        // Validate precision
        if (!Number.isInteger(this.options.precision) || this.options.precision < 1 || this.options.precision > 16) {
            throw new Error('Precision must be an integer between 1 and 16');
        }
        // Validate maxFileSize
        if (!Number.isInteger(this.options.maxFileSize) || this.options.maxFileSize <= 0) {
            throw new Error('Max file size must be a positive integer');
        }
        // Validate supportedFormats structure
        if (!this.options.supportedFormats || typeof this.options.supportedFormats !== 'object') {
            throw new Error('Invalid supported formats configuration');
        }
        for (const format in this.options.supportedFormats) {
            if (!Array.isArray(this.options.supportedFormats[format])) {
                throw new Error(`Invalid extensions array for format: ${ format }`);
            }
            for (const ext of this.options.supportedFormats[format]) {
                if (typeof ext !== 'string' || !ext.startsWith('.')) {
                    throw new Error(`Invalid file extension format: ${ ext }`);
                }
            }
        }
        // Validate boolean options
        if (typeof this.options.validateInput !== 'boolean') {
            throw new Error('validateInput must be a boolean');
        }
        if (typeof this.options.preserveMetadata !== 'boolean') {
            throw new Error('preserveMetadata must be a boolean');
        }
    }
    // Helper method to format IGES lines
    formatIGESLine(content, section, sequenceNumber) {
        // Pad content to 72 characters
        const paddedContent = content.padEnd(72, ' ');
        // Add section letter and sequence number
        return paddedContent + section + String(sequenceNumber).padStart(7, '0');
    }
    // Helper method to write NURBS curve
    writeNurbsCurve(curve, entityNumber, paramStart) {
        const directory = [];
        const parameter = [];
        // Entity type 126 (NURBS Curve)
        const paramData = [
            '126,',
            // Entity type
            curve.degree + ',',
            // Degree
            curve.controlPoints.length + ',',
            // Number of control points
            curve.knotVector.knots.length + ',',
            // Number of knots
            '0,',
            // Planar flag (0=non-planar)
            '0,',
            // Closed flag (0=open)
            '0,',
            // Rational flag (0=non-rational)
            '0,'
        ].join('');
        // Add knot vector
        const knotData = curve.knotVector.knots.map(k => k.toFixed(settings.precision)).join(',') + ',';
        // Add control points
        const controlPointData = curve.controlPoints.map(cp => {
            const [x, y, z] = cp.position();
            const w = cp.weight();
            return `${ x.toFixed(settings.precision) },${ y.toFixed(settings.precision) },${ z.toFixed(settings.precision) },${ w.toFixed(settings.precision) },`;
        }).join('');
        // Split parameter data into lines of appropriate length
        const combinedData = paramData + knotData + controlPointData;
        const parameterLines = this.splitParameterData(combinedData);
        parameter.push(...parameterLines);
        // Create directory entry
        directory.push(this.createDirectoryEntry(126, paramStart, entityNumber));
        return {
            directory,
            parameter
        };
    }
    // Helper method to write NURBS surface
    writeNurbsSurface(surface, entityNumber, paramStart) {
        const directory = [];
        const parameter = [];
        // Entity type 128 (NURBS Surface)
        const paramData = [
            '128,',
            // Entity type
            surface.degreeU + ',',
            // U degree
            surface.degreeV + ',',
            // V degree
            surface.controlPoints.length + ',',
            // Number of U control points
            surface.controlPoints[0].length + ',',
            // Number of V control points
            surface.knotVectorU.knots.length + ',',
            // Number of U knots
            surface.knotVectorV.knots.length + ',',
            // Number of V knots
            '0,',
            // U closed flag
            '0,'
        ].join('');
        // Add knot vectors
        const knotDataU = surface.knotVectorU.knots.map(k => k.toFixed(settings.precision)).join(',') + ',';
        const knotDataV = surface.knotVectorV.knots.map(k => k.toFixed(settings.precision)).join(',') + ',';
        // Add control points
        const controlPointData = surface.controlPoints.flat().map(cp => {
            const [x, y, z] = cp.position();
            const w = cp.weight();
            return `${ x.toFixed(settings.precision) },${ y.toFixed(settings.precision) },${ z.toFixed(settings.precision) },${ w.toFixed(settings.precision) },`;
        }).join('');
        // Split parameter data into lines of appropriate length
        const combinedData = paramData + knotDataU + knotDataV + controlPointData;
        const parameterLines = this.splitParameterData(combinedData);
        parameter.push(...parameterLines);
        // Create directory entry
        directory.push(this.createDirectoryEntry(128, paramStart, entityNumber));
        return {
            directory,
            parameter
        };
    }
    // Helper method to split parameter data into lines
    splitParameterData(data) {
        const lines = [];
        let remaining = data;
        const maxLength = 64;
        // Maximum characters per parameter line
        while (remaining.length > 0) {
            const line = remaining.slice(0, maxLength);
            remaining = remaining.slice(maxLength);
            lines.push(this.formatIGESLine(line, 'P', lines.length + 1));
        }
        return lines;
    }
    // Helper method to create directory entry
    createDirectoryEntry(entityType, paramStart, entityNumber) {
        // First directory entry line
        const line1 = [
            String(entityType).padStart(8),
            // Entity type
            String(paramStart).padStart(8),
            // Parameter pointer
            '0'.padStart(8),
            // Structure
            '0'.padStart(8),
            // Line Font Pattern
            '0'.padStart(8),
            // Level
            '0'.padStart(8),
            // View
            '0'.padStart(8),
            // Transform Matrix
            '0'.padStart(8),
            // Label Display
            '00000000'
        ].join('');
        // Second directory entry line
        const line2 = [
            String(entityType).padStart(8),
            // Entity type
            '0'.padStart(8),
            // Line Weight
            '0'.padStart(8),
            // Color
            '0'.padStart(8),
            // Parameter Line Count
            '0'.padStart(8),
            // Form Number
            '0'.padStart(8),
            // Reserved
            '0'.padStart(8),
            // Reserved
            '0'.padStart(8),
            // Entity Label
            '00000000'
        ].join('');
        return [
            this.formatIGESLine(line1, 'D', entityNumber * 2 - 1),
            this.formatIGESLine(line2, 'D', entityNumber * 2)
        ];
    }
    parseEntityId(line) {
        const match = line.match(/^#(\d+)/);
        if (!match) {
            throw new Error('Invalid entity ID format');
        }
        return parseInt(match[1]);
    }
    parseEntityType(line) {
        const match = line.match(/^#\d+\s*=\s*([A-Z_]+)/);
        if (!match) {
            throw new Error('Invalid entity type format');
        }
        return match[1];
    }
    parseHeaderAttribute(line) {
        const match = line.match(/'([^']+)'/);
        return match ? match[1] : '';
    }
    parseTimeStamp(line) {
        const match = line.match(/'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}'/);
        return match ? new Date(match[0].replace(/'/g, '')) : null;
    }
    processEntity(entity, entities) {
        const {id, type, data} = entity;
        // Parse parameters from entity data
        const params = this.parseParameters(data);
        // Store processed entity
        entities.set(id, {
            type,
            parameters: params,
            processed: false
        });
    }
    parseParameters(data) {
        // Remove entity header
        const paramString = data.replace(/^#\d+\s*=\s*[A-Z_]+\s*\(/, '').replace(/\);$/, '');
        // Split parameters and process each one
        return paramString.split(',').map(param => param.trim()).map(param => {
            if (param.startsWith('#')) {
                // Reference to another entity
                return {
                    type: 'reference',
                    value: parseInt(param.substring(1))
                };
            } else if (param.startsWith('(') && param.endsWith(')')) {
                // Parameter list
                return {
                    type: 'list',
                    value: this.parseParameters(param.substring(1, param.length - 1))
                };
            } else if (param.startsWith('\'') && param.endsWith('\'')) {
                // String value
                return {
                    type: 'string',
                    value: param.substring(1, param.length - 1)
                };
            } else if (!isNaN(param)) {
                // Numeric value
                return {
                    type: 'number',
                    value: parseFloat(param)
                };
            } else {
                // Enumeration or other value
                return {
                    type: 'enum',
                    value: param
                };
            }
        });
    }
    convertToGeometry(entity, entities) {
        switch (entity.type) {
        case 'B_SPLINE_CURVE_WITH_KNOTS':
            return this.convertBSplineCurve(entity, entities);
        case 'B_SPLINE_SURFACE_WITH_KNOTS':
            return this.convertBSplineSurface(entity, entities);
        case 'TRIMMED_CURVE':
            return this.convertTrimmedCurve(entity, entities);
        // Add other geometry types as needed
        default:
            return null;
        }
    }
    writeNurbsCurveSTEP(curve, startId, settings) {
        const entities = [];
        let currentId = startId;
        // Write control points as cartesian points
        const controlPointIds = curve.controlPoints.map(cp => {
            const [x, y, z] = cp.position();
            entities.push(`#${ currentId }=CARTESIAN_POINT('',` + `(${ x.toFixed(settings.precision) },` + `${ y.toFixed(settings.precision) },` + `${ z.toFixed(settings.precision) }));`);
            return currentId++;
        });
        // Write weights if curve is rational
        const weights = curve.controlPoints.map(cp => cp.weight());
        const isRational = weights.some(w => Math.abs(w - 1) > Number.EPSILON);
        let weightsRef = '';
        if (isRational) {
            const weightsStr = weights.map(w => w.toFixed(settings.precision)).join(',');
            entities.push(`#${ currentId }=B_SPLINE_CURVE_WITH_WEIGHTS(${ weightsStr });`);
            weightsRef = `#${ currentId }`;
            currentId++;
        }
        // Write knot multiplicities
        const knots = curve.knotVector.knots;
        const multiplicities = [];
        let currentKnot = knots[0];
        let currentMultiplicity = 1;
        for (let i = 1; i <= knots.length; i++) {
            if (i < knots.length && Math.abs(knots[i] - currentKnot) < Number.EPSILON) {
                currentMultiplicity++;
            } else {
                multiplicities.push(currentMultiplicity);
                if (i < knots.length) {
                    currentKnot = knots[i];
                    currentMultiplicity = 1;
                }
            }
        }
        // Write knots
        const uniqueKnots = [];
        for (let i = 0; i < knots.length; i++) {
            if (i === 0 || Math.abs(knots[i] - knots[i - 1]) >= Number.EPSILON) {
                uniqueKnots.push(knots[i]);
            }
        }
        // Write B-Spline curve entity
        const curveEntity = `#${ currentId }=B_SPLINE_CURVE_WITH_KNOTS('',` + `${ curve.degree },` + // Degree
        `(${ controlPointIds.join(',') }),` + // Control points
        (isRational ? `${ weightsRef },` : '') + // Weights reference if rational
        `B_SPLINE_CURVE_FORM(.UNSPECIFIED.),` + `.F.,` + // Not closed
        `.F.,` + // Not self intersecting
        `(${ multiplicities.join(',') }),` + // Knot multiplicities
        `(${ uniqueKnots.map(k => k.toFixed(settings.precision)).join(',') }),` + // Knot values
        `.UNSPECIFIED.);`;
        // Knot type
        entities.push(curveEntity);
        currentId++;
        return {
            data: entities,
            nextId: currentId
        };
    }
    writeNurbsSurfaceSTEP(surface, startId, settings) {
        const entities = [];
        let currentId = startId;
        // Write control points as cartesian points
        const controlPointIds = [];
        for (let i = 0; i < surface.controlPoints.length; i++) {
            const row = [];
            for (let j = 0; j < surface.controlPoints[i].length; j++) {
                const cp = surface.controlPoints[i][j];
                const [x, y, z] = cp.position();
                entities.push(`#${ currentId }=CARTESIAN_POINT('',` + `(${ x.toFixed(settings.precision) },` + `${ y.toFixed(settings.precision) },` + `${ z.toFixed(settings.precision) }));`);
                row.push(currentId++);
            }
            controlPointIds.push(row);
        }
        // Write weights if surface is rational
        const weights = surface.controlPoints.map(row => row.map(cp => cp.weight()));
        const isRational = weights.some(row => row.some(w => Math.abs(w - 1) > Number.EPSILON));
        let weightsRef = '';
        if (isRational) {
            const weightsStr = weights.map(row => `(${ row.map(w => w.toFixed(settings.precision)).join(',') })`).join(',');
            entities.push(`#${ currentId }=B_SPLINE_SURFACE_WITH_WEIGHTS(${ weightsStr });`);
            weightsRef = `#${ currentId }`;
            currentId++;
        }
        // Process U and V knot vectors
        const processKnotVector = knots => {
            const multiplicities = [];
            const uniqueKnots = [];
            let currentKnot = knots[0];
            let currentMultiplicity = 1;
            for (let i = 1; i <= knots.length; i++) {
                if (i < knots.length && Math.abs(knots[i] - currentKnot) < Number.EPSILON) {
                    currentMultiplicity++;
                } else {
                    multiplicities.push(currentMultiplicity);
                    uniqueKnots.push(currentKnot);
                    if (i < knots.length) {
                        currentKnot = knots[i];
                        currentMultiplicity = 1;
                    }
                }
            }
            return {
                multiplicities,
                uniqueKnots
            };
        };
        const uKnots = processKnotVector(surface.knotVectorU.knots);
        const vKnots = processKnotVector(surface.knotVectorV.knots);
        // Write B-Spline surface entity
        const surfaceEntity = `#${ currentId }=B_SPLINE_SURFACE_WITH_KNOTS('',` + `${ surface.degreeU },${ surface.degreeV },` + // Degrees
        `${ surface.controlPoints.length },${ surface.controlPoints[0].length },` + // Sizes
        `((${ controlPointIds.map(row => row.join(',')).join('),(') })),` + // Control points
        (isRational ? `${ weightsRef },` : '') + // Weights reference if rational
        `B_SPLINE_SURFACE_FORM(.UNSPECIFIED.),` + `.F.,.F.,.F.,.F.,` + // Form flags
        `(${ uKnots.multiplicities.join(',') }),` + // U multiplicities
        `(${ vKnots.multiplicities.join(',') }),` + // V multiplicities
        `(${ uKnots.uniqueKnots.map(k => k.toFixed(settings.precision)).join(',') }),` + // U knots
        `(${ vKnots.uniqueKnots.map(k => k.toFixed(settings.precision)).join(',') }),` + // V knots
        `.UNSPECIFIED.);`;
        // Surface type
        entities.push(surfaceEntity);
        currentId++;
        return {
            data: entities,
            nextId: currentId
        };
    }
    writeTrimmedSurfaceSTEP(surface, startId, settings) {
        const entities = [];
        let currentId = startId;
        // First write the base surface
        const baseSurfaceEntities = this.writeNurbsSurfaceSTEP(surface.baseSurface, currentId, settings);
        entities.push(...baseSurfaceEntities.data);
        const baseSurfaceId = baseSurfaceEntities.nextId - 1;
        currentId = baseSurfaceEntities.nextId;
        // Write trim curves
        const trimCurveIds = [];
        for (const trimLoop of surface.trimCurves) {
            const loopCurveIds = [];
            for (const curve of trimLoop.curves) {
                const curveEntities = this.writeNurbsCurveSTEP(curve, currentId, settings);
                entities.push(...curveEntities.data);
                loopCurveIds.push(curveEntities.nextId - 1);
                currentId = curveEntities.nextId;
            }
            trimCurveIds.push({
                curves: loopCurveIds,
                isOuter: trimLoop.isOuter
            });
        }
        // Write trimmed surface entity
        const trimmedEntity = `#${ currentId }=TRIMMED_SURFACE('',` + `#${ baseSurfaceId },` + // Base surface reference
        `(${ trimCurveIds.map(loop => `(${ loop.curves.join(',') })`).join(',') }),` + // Trim curves
        `${ trimCurveIds.some(loop => loop.isOuter) },` + // Has outer boundary
        `.T.);`;
        // Sense agreement flag
        entities.push(trimmedEntity);
        currentId++;
        return {
            data: entities,
            nextId: currentId
        };
    }
    parseNurbsCurve(paramData, globalData) {
        try {
            // Validate input
            if (!Array.isArray(paramData)) {
                throw new Error('Parameter data must be an array');
            }
            if (!globalData || typeof globalData !== 'object') {
                throw new Error('Invalid global data');
            }
            // Extract parameters
            const degree = parseInt(paramData[1]);
            const numPoles = parseInt(paramData[2]);
            const numKnots = parseInt(paramData[3]);
            const isRational = parseInt(paramData[6]) !== 0;
            // Validate degree and counts
            if (!Number.isInteger(degree) || degree < 1) {
                throw new Error('Invalid curve degree');
            }
            if (!Number.isInteger(numPoles) || numPoles < degree + 1) {
                throw new Error('Invalid number of control points');
            }
            if (!Number.isInteger(numKnots) || numKnots !== numPoles + degree + 1) {
                throw new Error('Invalid number of knots');
            }
            // Extract knot vector (starting at index 8)
            const knots = [];
            let currentIndex = 8;
            for (let i = 0; i < numKnots; i++) {
                const knotValue = parseFloat(paramData[currentIndex++]);
                if (!Number.isFinite(knotValue)) {
                    throw new Error(`Invalid knot value at index ${ i }`);
                }
                knots.push(knotValue);
            }
            // Extract weights if curve is rational
            const weights = [];
            if (isRational) {
                for (let i = 0; i < numPoles; i++) {
                    const weight = parseFloat(paramData[currentIndex++]);
                    if (!Number.isFinite(weight) || weight <= 0) {
                        throw new Error(`Invalid weight value at index ${ i }`);
                    }
                    weights.push(weight);
                }
            }
            // Extract control points
            const controlPoints = [];
            for (let i = 0; i < numPoles; i++) {
                const x = parseFloat(paramData[currentIndex++]) * globalData.modelSpaceScale;
                const y = parseFloat(paramData[currentIndex++]) * globalData.modelSpaceScale;
                const z = parseFloat(paramData[currentIndex++]) * globalData.modelSpaceScale;
                if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
                    throw new Error(`Invalid control point coordinates at index ${ i }`);
                }
                const weight = isRational ? weights[i] : 1;
                controlPoints.push(new ControlPoint(x, y, z, weight));
            }
            // Create knot vector
            const knotVector = new KnotVector(knots);
            // Validate control points and knot vector
            if (controlPoints.length !== numPoles) {
                throw new Error('Incorrect number of control points parsed');
            }
            if (knotVector.knots.length !== numKnots) {
                throw new Error('Incorrect number of knots parsed');
            }
            // Create and return new NURBS curve
            return new NurbsCurve(controlPoints, knotVector, degree);
        } catch (error) {
            throw new Error(`Failed to parse NURBS curve: ${ error.message }`);
        }
    }
    parseNurbsSurface(paramData, globalData) {
        try {
            // Validate input
            if (!Array.isArray(paramData)) {
                throw new Error('Parameter data must be an array');
            }
            if (!globalData || typeof globalData !== 'object') {
                throw new Error('Invalid global data');
            }
            // Extract parameters
            const degreeU = parseInt(paramData[1]);
            const degreeV = parseInt(paramData[2]);
            const numPolesU = parseInt(paramData[3]);
            const numPolesV = parseInt(paramData[4]);
            const numKnotsU = parseInt(paramData[5]);
            const numKnotsV = parseInt(paramData[6]);
            const isRational = parseInt(paramData[7]) !== 0;
            // Validate degrees and counts
            if (!Number.isInteger(degreeU) || !Number.isInteger(degreeV) || degreeU < 1 || degreeV < 1) {
                throw new Error('Invalid surface degrees');
            }
            if (!Number.isInteger(numPolesU) || !Number.isInteger(numPolesV) || numPolesU < degreeU + 1 || numPolesV < degreeV + 1) {
                throw new Error('Invalid number of control points');
            }
            if (!Number.isInteger(numKnotsU) || !Number.isInteger(numKnotsV) || numKnotsU !== numPolesU + degreeU + 1 || numKnotsV !== numPolesV + degreeV + 1) {
                throw new Error('Invalid number of knots');
            }
            // Extract U knot vector (starting at index 8)
            const knotsU = [];
            let currentIndex = 8;
            for (let i = 0; i < numKnotsU; i++) {
                const knotValue = parseFloat(paramData[currentIndex++]);
                if (!Number.isFinite(knotValue)) {
                    throw new Error(`Invalid U knot value at index ${ i }`);
                }
                knotsU.push(knotValue);
            }
            // Extract V knot vector
            const knotsV = [];
            for (let i = 0; i < numKnotsV; i++) {
                const knotValue = parseFloat(paramData[currentIndex++]);
                if (!Number.isFinite(knotValue)) {
                    throw new Error(`Invalid V knot value at index ${ i }`);
                }
                knotsV.push(knotValue);
            }
            // Extract weights if surface is rational
            const weights = [];
            if (isRational) {
                for (let i = 0; i < numPolesU; i++) {
                    const rowWeights = [];
                    for (let j = 0; j < numPolesV; j++) {
                        const weight = parseFloat(paramData[currentIndex++]);
                        if (!Number.isFinite(weight) || weight <= 0) {
                            throw new Error(`Invalid weight value at index (${ i },${ j })`);
                        }
                        rowWeights.push(weight);
                    }
                    weights.push(rowWeights);
                }
            }
            // Extract control points
            const controlPoints = Array(numPolesU).fill().map(() => Array(numPolesV));
            for (let i = 0; i < numPolesU; i++) {
                for (let j = 0; j < numPolesV; j++) {
                    const x = parseFloat(paramData[currentIndex++]) * globalData.modelSpaceScale;
                    const y = parseFloat(paramData[currentIndex++]) * globalData.modelSpaceScale;
                    const z = parseFloat(paramData[currentIndex++]) * globalData.modelSpaceScale;
                    if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
                        throw new Error(`Invalid control point coordinates at index (${ i },${ j })`);
                    }
                    const weight = isRational ? weights[i][j] : 1;
                    controlPoints[i][j] = new ControlPoint(x, y, z, weight);
                }
            }
            // Create knot vectors
            const knotVectorU = new KnotVector(knotsU);
            const knotVectorV = new KnotVector(knotsV);
            // Create and validate surface dimensions
            if (controlPoints.length !== numPolesU || !controlPoints.every(row => row.length === numPolesV)) {
                throw new Error('Incorrect control point grid dimensions');
            }
            // Create and return new NURBS surface
            const surface = new NurbsSurface();
            surface.controlPoints = controlPoints;
            surface.knotVectorU = knotVectorU;
            surface.knotVectorV = knotVectorV;
            surface.degreeU = degreeU;
            surface.degreeV = degreeV;
            return surface;
        } catch (error) {
            throw new Error(`Failed to parse NURBS surface: ${ error.message }`);
        }
    }
    parseTrimmedSurface(paramData, globalData, entities, parameters) {
        try {
            // Validate inputs
            if (!Array.isArray(paramData)) {
                throw new Error('Parameter data must be an array');
            }
            if (!globalData || typeof globalData !== 'object') {
                throw new Error('Invalid global data');
            }
            // Extract base surface reference (first parameter)
            const baseSurfaceIndex = parseInt(paramData[1]);
            if (!Number.isInteger(baseSurfaceIndex) || baseSurfaceIndex <= 0) {
                throw new Error('Invalid base surface reference');
            }
            // Get base surface data
            const baseSurfaceParams = parameters.get(baseSurfaceIndex);
            if (!baseSurfaceParams) {
                throw new Error('Base surface parameters not found');
            }
            // Parse base surface
            const baseSurface = this.parseNurbsSurface(baseSurfaceParams, globalData);
            if (!baseSurface) {
                throw new Error('Failed to parse base surface');
            }
            // Create trimmed surface instance
            const trimmedSurface = new TrimmedSurface();
            trimmedSurface.baseSurface = baseSurface;
            // Parse outer boundary curves (third parameter)
            const outerBoundaryIndex = parseInt(paramData[3]);
            if (Number.isInteger(outerBoundaryIndex) && outerBoundaryIndex > 0) {
                const outerBoundaryParams = parameters.get(outerBoundaryIndex);
                if (outerBoundaryParams) {
                    const outerCurves = this.parseTrimCurves(outerBoundaryParams, globalData);
                    if (outerCurves && outerCurves.length > 0) {
                        trimmedSurface.setTrimCurves(outerCurves, true);
                    }
                }
            }
            // true for outer boundary
            // Parse inner boundary curves (fifth parameter)
            const innerBoundaryCount = parseInt(paramData[5]);
            if (Number.isInteger(innerBoundaryCount) && innerBoundaryCount > 0) {
                for (let i = 0; i < innerBoundaryCount; i++) {
                    const innerBoundaryIndex = parseInt(paramData[6 + i]);
                    if (Number.isInteger(innerBoundaryIndex) && innerBoundaryIndex > 0) {
                        const innerBoundaryParams = parameters.get(innerBoundaryIndex);
                        if (innerBoundaryParams) {
                            const innerCurves = this.parseTrimCurves(innerBoundaryParams, globalData);
                            if (innerCurves && innerCurves.length > 0) {
                                trimmedSurface.setTrimCurves(innerCurves, false);
                            }
                        }
                    }
                }
            }
            // false for inner boundary
            // Validate the trimmed surface
            if (!trimmedSurface.baseSurface) {
                throw new Error('Invalid trimmed surface: missing base surface');
            }
            if (trimmedSurface.trimCurves.length === 0) {
                throw new Error('Invalid trimmed surface: no trim curves defined');
            }
            // Return the parsed trimmed surface
            return trimmedSurface;
        } catch (error) {
            throw new Error(`Failed to parse trimmed surface: ${ error.message }`);
        }
    }
    convertBSplineCurve(entity, entities) {
        try {
            // Validate entity type
            if (entity.type !== 'B_SPLINE_CURVE_WITH_KNOTS') {
                throw new Error('Invalid entity type for B-spline curve conversion');
            }
            // Extract parameters from entity
            const degree = entity.parameters[0].value;
            const controlPointRefs = entity.parameters[1].value;
            const knotMultiplicities = entity.parameters[6].value;
            const knots = entity.parameters[7].value;
            let weights = null;
            // Check if the curve is rational
            if (entity.parameters[2] && entity.parameters[2].type === 'reference') {
                const weightsEntity = entities.get(entity.parameters[2].value);
                if (weightsEntity && weightsEntity.parameters) {
                    weights = weightsEntity.parameters.map(w => w.value);
                }
            }
            // Convert control point references to actual points
            const controlPoints = [];
            for (const ref of controlPointRefs) {
                if (ref.type !== 'reference') {
                    throw new Error('Invalid control point reference');
                }
                const pointEntity = entities.get(ref.value);
                if (!pointEntity || !pointEntity.parameters) {
                    throw new Error('Invalid control point entity');
                }
                const coordinates = pointEntity.parameters.map(coord => coord.value);
                if (coordinates.length !== 3) {
                    throw new Error('Invalid control point coordinates');
                }
                // Create control point with weight if curve is rational
                const weight = weights ? weights[controlPoints.length] : 1;
                controlPoints.push(new ControlPoint(coordinates[0], coordinates[1], coordinates[2], weight));
            }
            // Build knot vector from multiplicities and knot values
            const expandedKnots = [];
            for (let i = 0; i < knots.length; i++) {
                const multiplicity = knotMultiplicities[i].value;
                const knotValue = knots[i].value;
                for (let j = 0; j < multiplicity; j++) {
                    expandedKnots.push(knotValue);
                }
            }
            // Validate knot vector length
            const expectedKnotLength = controlPoints.length + degree + 1;
            if (expandedKnots.length !== expectedKnotLength) {
                throw new Error(`Invalid knot vector length. Expected ${ expectedKnotLength }, got ${ expandedKnots.length }`);
            }
            // Create knot vector
            const knotVector = new KnotVector(expandedKnots);
            // Create and return NURBS curve
            return new NurbsCurve(controlPoints, knotVector, degree);
        } catch (error) {
            throw new Error(`Failed to convert B-spline curve: ${ error.message }`);
        }
    }
    convertBSplineSurface(entity, entities) {
        try {
            // Validate entity type
            if (entity.type !== 'B_SPLINE_SURFACE_WITH_KNOTS') {
                throw new Error('Invalid entity type for B-spline surface conversion');
            }
            // Extract surface parameters
            const degreeU = entity.parameters[0].value;
            const degreeV = entity.parameters[1].value;
            const numControlPointsU = entity.parameters[2].value;
            const numControlPointsV = entity.parameters[3].value;
            const controlPointRefs = entity.parameters[4].value;
            const knotMultiplicitiesU = entity.parameters[7].value;
            const knotMultiplicitiesV = entity.parameters[8].value;
            const knotsU = entity.parameters[9].value;
            const knotsV = entity.parameters[10].value;
            // Check for rational surface (weights)
            let weights = null;
            if (entity.parameters[5] && entity.parameters[5].type === 'reference') {
                const weightsEntity = entities.get(entity.parameters[5].value);
                if (weightsEntity && weightsEntity.parameters) {
                    weights = weightsEntity.parameters.map(row => row.value.map(w => w.value));
                }
            }
            // Convert control point references to actual points
            const controlPoints = Array(numControlPointsU).fill().map(() => Array(numControlPointsV).fill(null));
            // Process control point grid
            for (let i = 0; i < numControlPointsU; i++) {
                for (let j = 0; j < numControlPointsV; j++) {
                    const ref = controlPointRefs[i][j];
                    if (ref.type !== 'reference') {
                        throw new Error(`Invalid control point reference at (${ i }, ${ j })`);
                    }
                    const pointEntity = entities.get(ref.value);
                    if (!pointEntity || !pointEntity.parameters) {
                        throw new Error(`Invalid control point entity at (${ i }, ${ j })`);
                    }
                    const coordinates = pointEntity.parameters.map(coord => coord.value);
                    if (coordinates.length !== 3) {
                        throw new Error(`Invalid control point coordinates at (${ i }, ${ j })`);
                    }
                    // Create control point with weight if surface is rational
                    const weight = weights ? weights[i][j] : 1;
                    controlPoints[i][j] = new ControlPoint(coordinates[0], coordinates[1], coordinates[2], weight);
                }
            }
            // Build knot vectors from multiplicities and knot values
            const expandedKnotsU = [];
            for (let i = 0; i < knotsU.length; i++) {
                const multiplicity = knotMultiplicitiesU[i].value;
                const knotValue = knotsU[i].value;
                for (let j = 0; j < multiplicity; j++) {
                    expandedKnotsU.push(knotValue);
                }
            }
            const expandedKnotsV = [];
            for (let i = 0; i < knotsV.length; i++) {
                const multiplicity = knotMultiplicitiesV[i].value;
                const knotValue = knotsV[i].value;
                for (let j = 0; j < multiplicity; j++) {
                    expandedKnotsV.push(knotValue);
                }
            }
            // Validate knot vector lengths
            const expectedKnotLengthU = numControlPointsU + degreeU + 1;
            const expectedKnotLengthV = numControlPointsV + degreeV + 1;
            if (expandedKnotsU.length !== expectedKnotLengthU) {
                throw new Error(`Invalid U knot vector length. Expected ${ expectedKnotLengthU }, got ${ expandedKnotsU.length }`);
            }
            if (expandedKnotsV.length !== expectedKnotLengthV) {
                throw new Error(`Invalid V knot vector length. Expected ${ expectedKnotLengthV }, got ${ expandedKnotsV.length }`);
            }
            // Create knot vectors
            const knotVectorU = new KnotVector(expandedKnotsU);
            const knotVectorV = new KnotVector(expandedKnotsV);
            // Create and return NURBS surface
            const surface = new NurbsSurface();
            surface.controlPoints = controlPoints;
            surface.knotVectorU = knotVectorU;
            surface.knotVectorV = knotVectorV;
            surface.degreeU = degreeU;
            surface.degreeV = degreeV;
            return surface;
        } catch (error) {
            throw new Error(`Failed to convert B-spline surface: ${ error.message }`);
        }
    }
    convertTrimmedCurve(entity, entities) {
        try {
            // Validate entity type
            if (entity.type !== 'TRIMMED_CURVE') {
                throw new Error('Invalid entity type for trimmed curve conversion');
            }
            // Extract base curve reference and trim parameters
            const baseCurveRef = entity.parameters[0];
            if (!baseCurveRef || baseCurveRef.type !== 'reference') {
                throw new Error('Invalid base curve reference');
            }
            // Get base curve entity
            const baseCurveEntity = entities.get(baseCurveRef.value);
            if (!baseCurveEntity) {
                throw new Error('Base curve entity not found');
            }
            // Convert base curve
            let baseCurve;
            if (baseCurveEntity.type === 'B_SPLINE_CURVE_WITH_KNOTS') {
                baseCurve = this.convertBSplineCurve(baseCurveEntity, entities);
            } else {
                throw new Error('Unsupported base curve type');
            }
            // Extract trim parameters
            const trimParams = entity.parameters[1];
            if (!Array.isArray(trimParams.value) || trimParams.value.length !== 2) {
                throw new Error('Invalid trim parameters');
            }
            // Get start and end parameters
            const startParam = this.extractTrimParameter(trimParams.value[0], entities);
            const endParam = this.extractTrimParameter(trimParams.value[1], entities);
            if (typeof startParam !== 'number' || typeof endParam !== 'number') {
                throw new Error('Invalid trim parameter values');
            }
            // Check if parameters are within base curve domain
            if (startParam < baseCurve.domain.min || startParam > baseCurve.domain.max || endParam < baseCurve.domain.min || endParam > baseCurve.domain.max) {
                throw new Error('Trim parameters outside curve domain');
            }
            // Split the curve at trim points
            const midCurves = baseCurve.split(startParam);
            if (!midCurves || midCurves.length !== 2) {
                throw new Error('Failed to split curve at start parameter');
            }
            const trimmedPortion = midCurves[1].split(endParam)[0];
            if (!trimmedPortion) {
                throw new Error('Failed to split curve at end parameter');
            }
            // Check sense agreement (reverse if necessary)
            const senseAgreement = entity.parameters[2].value;
            if (typeof senseAgreement === 'boolean' && !senseAgreement) {
                return trimmedPortion.reverse();
            }
            return trimmedPortion;
        } catch (error) {
            throw new Error(`Failed to convert trimmed curve: ${ error.message }`);
        }
    }
    extractTrimParameter(paramEntity, entities) {
        // Handle different types of trim parameters
        if (paramEntity.type === 'number') {
            return paramEntity.value;
        }
        if (paramEntity.type === 'reference') {
            const refEntity = entities.get(paramEntity.value);
            if (!refEntity) {
                throw new Error('Referenced parameter entity not found');
            }
            // Handle CARTESIAN_POINT
            if (refEntity.type === 'CARTESIAN_POINT') {
                // Convert point to parameter using base curve
                // This would require additional implementation for point inversion
                throw new Error('Point-based trim parameters not yet supported');
            }
            // Handle PARAMETER_VALUE
            if (refEntity.type === 'PARAMETER_VALUE') {
                return refEntity.parameters[0].value;
            }
        }
        throw new Error('Unsupported trim parameter type');
    }
    // Helper method to parse trim curves
    parseTrimCurves(paramData, globalData) {
        try {
            const curves = [];
            let currentIndex = 1;
            // Skip first parameter
            while (currentIndex < paramData.length) {
                const curveType = parseInt(paramData[currentIndex++]);
                switch (curveType) {
                case 126:
                    // NURBS curve
                    const curveParams = [];
                    // Collect parameters for this curve
                    while (currentIndex < paramData.length && !Number.isInteger(parseInt(paramData[currentIndex]))) {
                        curveParams.push(paramData[currentIndex++]);
                    }
                    if (curveParams.length > 0) {
                        const curve = this.parseNurbsCurve(curveParams, globalData);
                        if (curve) {
                            curves.push(curve);
                        }
                    }
                    break;
                // Add other curve types as needed
                default:
                    // Skip unknown curve type
                    while (currentIndex < paramData.length && !Number.isInteger(parseInt(paramData[currentIndex]))) {
                        currentIndex++;
                    }
                    break;
                }
            }
            // Validate that curves form a closed loop
            for (let i = 0; i < curves.length; i++) {
                const currentCurve = curves[i];
                const nextCurve = curves[(i + 1) % curves.length];
                const endPoint = currentCurve.evaluate(currentCurve.domain.max);
                const startPoint = nextCurve.evaluate(nextCurve.domain.min);
                if (endPoint.subtract(startPoint).length() > this.tolerance) {
                    throw new Error('Trim curves do not form a closed loop');
                }
            }
            return curves;
        } catch (error) {
            throw new Error(`Failed to parse trim curves: ${ error.message }`);
        }
    }
}
class GeometryContainer {
    constructor() {
        this.entities = new Map();
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
            throw new Error(`Geometry with ID "${ id }" already exists`);
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
            throw new Error(`Failed to add geometry: ${ error.message }`);
        }
    }
    remove(id) {
        // Input validation
        if (typeof id !== 'string' || id.trim().length === 0) {
            throw new Error('ID must be a non-empty string');
        }
        // Check if geometry exists
        if (!this.geometries.has(id)) {
            throw new Error(`Geometry with ID "${ id }" does not exist`);
        }
        try {
            // Get the geometry before removal for return value
            const geometry = this.geometries.get(id);
            // Remove the geometry
            const success = this.geometries.delete(id);
            if (!success) {
                throw new Error(`Failed to remove geometry with ID "${ id }"`);
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
            throw new Error(`Error removing geometry: ${ error.message }`);
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
                        throw new Error(`Invalid search criterion: ${ key }`);
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
            throw new Error(`Search failed: ${ error.message }`);
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
            throw new Error(`Geometry with ID "${ id }" does not exist`);
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
            throw new Error(`Transformation failed: ${ error.message }`);
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
                    validationResults.errors.push(`Invalid ID format for geometry: ${ id }`);
                    validationResults.isValid = false;
                    continue;
                }
                // Check geometry object structure
                if (!geometry || !geometry.type || !geometry.data) {
                    validationResults.errors.push(`Invalid geometry structure for ID: ${ id }`);
                    validationResults.isValid = false;
                    continue;
                }
                // Check geometry type
                if (![
                        'NurbsCurve',
                        'NurbsSurface',
                        'TrimmedSurface'
                    ].includes(geometry.type)) {
                    validationResults.errors.push(`Invalid geometry type for ID: ${ id }`);
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
                    validationResults.errors.push(`Invalid dateAdded for ID: ${ id }`);
                    validationResults.isValid = false;
                }
                if (!(geometry.lastModified instanceof Date) || isNaN(geometry.lastModified.getTime())) {
                    validationResults.errors.push(`Invalid lastModified for ID: ${ id }`);
                    validationResults.isValid = false;
                }
                // Check for future dates
                const now = new Date();
                if (geometry.dateAdded > now || geometry.lastModified > now) {
                    validationResults.warnings.push(`Future timestamps detected for ID: ${ id }`);
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
            validationResults.errors.push(`Validation failed: ${ error.message }`);
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
            results.errors.push(`Invalid NURBS curve instance for ID: ${ id }`);
            results.isValid = false;
            return;
        }
        // Validate control points
        if (!Array.isArray(curve.controlPoints) || curve.controlPoints.length < curve.degree + 1) {
            results.errors.push(`Invalid control points for curve ID: ${ id }`);
            results.isValid = false;
        }
        // Validate knot vector
        if (!(curve.knotVector instanceof KnotVector)) {
            results.errors.push(`Invalid knot vector for curve ID: ${ id }`);
            results.isValid = false;
        }
        // Validate degree
        if (!Number.isInteger(curve.degree) || curve.degree < 1) {
            results.errors.push(`Invalid degree for curve ID: ${ id }`);
            results.isValid = false;
        }
        // Check for degenerate curves
        try {
            const start = curve.evaluate(curve.domain.min);
            const end = curve.evaluate(curve.domain.max);
            if (start.subtract(end).length() < Number.EPSILON) {
                results.warnings.push(`Potentially degenerate curve detected for ID: ${ id }`);
            }
        } catch (error) {
            results.errors.push(`Curve evaluation failed for ID: ${ id }: ${ error.message }`);
            results.isValid = false;
        }
    }
    validateNurbsSurface(id, surface, results) {
        if (!(surface instanceof NurbsSurface)) {
            results.errors.push(`Invalid NURBS surface instance for ID: ${ id }`);
            results.isValid = false;
            return;
        }
        // Validate control point network
        if (!Array.isArray(surface.controlPoints) || surface.controlPoints.length === 0) {
            results.errors.push(`Invalid control point network for surface ID: ${ id }`);
            results.isValid = false;
        }
        // Check control point network regularity
        const width = surface.controlPoints[0].length;
        if (!surface.controlPoints.every(row => Array.isArray(row) && row.length === width)) {
            results.errors.push(`Irregular control point network for surface ID: ${ id }`);
            results.isValid = false;
        }
        // Validate knot vectors
        if (!(surface.knotVectorU instanceof KnotVector) || !(surface.knotVectorV instanceof KnotVector)) {
            results.errors.push(`Invalid knot vectors for surface ID: ${ id }`);
            results.isValid = false;
        }
        // Validate degrees
        if (!Number.isInteger(surface.degreeU) || !Number.isInteger(surface.degreeV) || surface.degreeU < 1 || surface.degreeV < 1) {
            results.errors.push(`Invalid degrees for surface ID: ${ id }`);
            results.isValid = false;
        }
        // Check for degenerate surfaces
        try {
            const area = this.calculateApproximateArea(surface);
            if (area < Number.EPSILON) {
                results.warnings.push(`Potentially degenerate surface detected for ID: ${ id }`);
            }
        } catch (error) {
            results.errors.push(`Surface evaluation failed for ID: ${ id }: ${ error.message }`);
            results.isValid = false;
        }
    }
    validateTrimmedSurface(id, surface, results) {
        if (!(surface instanceof TrimmedSurface)) {
            results.errors.push(`Invalid trimmed surface instance for ID: ${ id }`);
            results.isValid = false;
            return;
        }
        // Validate base surface
        if (!(surface.baseSurface instanceof NurbsSurface)) {
            results.errors.push(`Invalid base surface for trimmed surface ID: ${ id }`);
            results.isValid = false;
        }
        // Validate trim curves
        if (!Array.isArray(surface.trimCurves)) {
            results.errors.push(`Invalid trim curves for surface ID: ${ id }`);
            results.isValid = false;
            return;
        }
        // Check trim curve validity
        for (let i = 0; i < surface.trimCurves.length; i++) {
            const trimLoop = surface.trimCurves[i];
            if (!Array.isArray(trimLoop.curves) || !trimLoop.hasOwnProperty('isOuter')) {
                results.errors.push(`Invalid trim loop structure at index ${ i } for surface ID: ${ id }`);
                results.isValid = false;
                continue;
            }
            // Verify trim curves form closed loops
            const curves = trimLoop.curves;
            for (let j = 0; j < curves.length; j++) {
                if (!(curves[j] instanceof NurbsCurve)) {
                    results.errors.push(`Invalid trim curve at index ${ j } in loop ${ i } for surface ID: ${ id }`);
                    results.isValid = false;
                }
            }
        }
        // Check for valid outer boundary
        const hasOuterBoundary = surface.trimCurves.some(loop => loop.isOuter);
        if (!hasOuterBoundary && surface.trimCurves.length > 0) {
            results.errors.push(`Missing outer boundary for trimmed surface ID: ${ id }`);
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
    update(id, geometry) {
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
        // Check if ID exists
        if (!this.entities.has(id)) {
            throw new Error(`Geometry with ID "${ id }" does not exist`);
        }
        // Create a deep copy of the geometry to prevent external modifications
        let geometryCopy;
        if (geometry instanceof NurbsCurve) {
            geometryCopy = new NurbsCurve([...geometry.controlPoints], geometry.knotVector.clone(), geometry.degree);
        } else if (geometry instanceof NurbsSurface) {
            geometryCopy = new NurbsSurface();
            geometryCopy.controlPoints = geometry.controlPoints.map(row => [...row]);
            geometryCopy.knotVectorU = geometry.knotVectorU.clone();
            geometryCopy.knotVectorV = geometry.knotVectorV.clone();
            geometryCopy.degreeU = geometry.degreeU;
            geometryCopy.degreeV = geometry.degreeV;
        } else if (geometry instanceof TrimmedSurface) {
            geometryCopy = new TrimmedSurface();
            geometryCopy.baseSurface = geometry.baseSurface.clone();
            geometryCopy.trimCurves = geometry.trimCurves.map(loop => ({
                curves: [...loop.curves],
                isOuter: loop.isOuter
            }));
        }
        // Update the geometry in the container
        this.entities.set(id, {
            type: geometry.constructor.name,
            data: geometryCopy,
            dateAdded: this.entities.get(id).dateAdded,
            // Preserve original dateAdded
            lastModified: new Date()
        });
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
                    const pointKey = `${ point.x.toFixed(6) },${ point.y.toFixed(6) },${ point.z.toFixed(6) }`;
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
                        const pointKey = `${ p1.x.toFixed(6) },${ p1.y.toFixed(6) },${ p1.z.toFixed(6) }`;
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
                    const pointKey = `${ point.x.toFixed(6) },${ point.y.toFixed(6) },${ point.z.toFixed(6) }`;
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
            throw new Error(`Refinement iteration failed: ${ error.message }`);
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
    refineIntersection(curve1, curve2, t1, t2, tolerance) {
        // Input validation
        if (!(curve1 instanceof NurbsCurve) || !(curve2 instanceof NurbsCurve)) {
            throw new Error('Both inputs must be NURBS curves');
        }
        if (typeof t1 !== 'number' || typeof t2 !== 'number' || !Number.isFinite(t1) || !Number.isFinite(t2)) {
            throw new Error('Parameters t1 and t2 must be finite numbers');
        }
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        const maxIterations = 20;
        // Maximum number of Newton iterations
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
                // Calculate normal equations components
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
            // If no convergence after max iterations, return best result
            return {
                t1: currentT1,
                t2: currentT2,
                point: curve1.evaluate(currentT1),
                distance: curve1.evaluate(currentT1).subtract(curve2.evaluate(currentT2)).length()
            };
        } catch (error) {
            throw new Error(`Intersection refinement failed: ${ error.message }`);
        }
    }
    refineSurfaceCurveIntersection(surface, curve, u, v, t, tolerance = this.tolerance) {
        // Input validation
        if (!(surface instanceof NurbsSurface)) {
            throw new Error('First parameter must be a NURBS surface');
        }
        if (!(curve instanceof NurbsCurve)) {
            throw new Error('Second parameter must be a NURBS curve');
        }
        if (typeof u !== 'number' || typeof v !== 'number' || typeof t !== 'number') {
            throw new Error('Parameters u, v, and t must be numbers');
        }
        if (!Number.isFinite(u) || !Number.isFinite(v) || !Number.isFinite(t)) {
            throw new Error('Parameters u, v, and t must be finite numbers');
        }
        const maxIterations = 20;
        let currentU = u;
        let currentV = v;
        let currentT = t;
        try {
            // Newton-Raphson iteration
            for (let iteration = 0; iteration < maxIterations; iteration++) {
                // Get points on surface and curve
                const surfacePoint = surface.evaluate(currentU, currentV);
                const curvePoint = curve.evaluate(currentT);
                // Calculate distance vector
                const distance = curvePoint.subtract(surfacePoint);
                // Check if points are close enough
                if (distance.length() < tolerance) {
                    return {
                        u: currentU,
                        v: currentV,
                        t: currentT,
                        point: surfacePoint,
                        distance: distance.length()
                    };
                }
                // Calculate derivatives
                const surfaceU = surface.derivative(currentU, currentV, 1, 0);
                const surfaceV = surface.derivative(currentU, currentV, 0, 1);
                const curveDerivative = curve.derivative(currentT, 1);
                // Set up system of equations for Newton iteration
                // [Su.x   Sv.x   -C'.x] [du]   [dx]
                // [Su.y   Sv.y   -C'.y] [dv] = [dy]
                // [Su.z   Sv.z   -C'.z] [dt]   [dz]
                const A = [
                    [
                        surfaceU.x,
                        surfaceV.x,
                        -curveDerivative.x
                    ],
                    [
                        surfaceU.y,
                        surfaceV.y,
                        -curveDerivative.y
                    ],
                    [
                        surfaceU.z,
                        surfaceV.z,
                        -curveDerivative.z
                    ]
                ];
                const b = [
                    -distance.x,
                    -distance.y,
                    -distance.z
                ];
                // Solve system using LU decomposition
                const solution = this.solveLU(A, b);
                if (!solution) {
                    throw new Error('Failed to solve linear system');
                }
                const [du, dv, dt] = solution;
                // Apply damping to ensure convergence
                const damping = 0.5;
                const newU = currentU + damping * du;
                const newV = currentV + damping * dv;
                const newT = currentT + damping * dt;
                // Check for convergence
                if (Math.abs(newU - currentU) < tolerance && Math.abs(newV - currentV) < tolerance && Math.abs(newT - currentT) < tolerance) {
                    return {
                        u: newU,
                        v: newV,
                        t: newT,
                        point: surface.evaluate(newU, newV),
                        distance: surface.evaluate(newU, newV).subtract(curve.evaluate(newT)).length()
                    };
                }
                // Update current parameters
                currentU = Math.max(surface.knotVectorU.domain.min, Math.min(surface.knotVectorU.domain.max, newU));
                currentV = Math.max(surface.knotVectorV.domain.min, Math.min(surface.knotVectorV.domain.max, newV));
                currentT = Math.max(curve.domain.min, Math.min(curve.domain.max, newT));
            }
            // If no convergence after max iterations, return best result
            return {
                u: currentU,
                v: currentV,
                t: currentT,
                point: surface.evaluate(currentU, currentV),
                distance: surface.evaluate(currentU, currentV).subtract(curve.evaluate(currentT)).length()
            };
        } catch (error) {
            throw new Error(`Surface-curve intersection refinement failed: ${ error.message }`);
        }
    }
    calculate(curve1, curve2) {
        // Input validation
        if (!(curve1 instanceof NurbsCurve)) {
            throw new Error('First parameter must be a NURBS curve');
        }
        if (!(curve2 instanceof NurbsCurve)) {
            throw new Error('Second parameter must be a NURBS curve');
        }
        // Initialize result object
        const result = new IntersectionResult();
        result.points = [];
        result.parameters = [];
        // Check bounding boxes first
        const box1 = this.computeBoundingBox(curve1);
        const box2 = this.computeBoundingBox(curve2);
        if (!this.boundingBoxesIntersect(box1, box2)) {
            return result;
        }
        // Initialize subdivision intervals for both curves
        const intervals1 = [{
                t1: curve1.domain.min,
                t2: curve1.domain.max
            }];
        const intervals2 = [{
                t1: curve2.domain.min,
                t2: curve2.domain.max
            }];
        // Store candidate intersection points
        const candidatePoints = new Set();
        // Process subdivided intervals
        while (intervals1.length > 0 && intervals2.length > 0) {
            const i1 = intervals1.pop();
            const i2 = intervals2.pop();
            // Compute bounding boxes for current intervals
            const box1Current = this.computeIntervalBoundingBox(curve1, i1.t1, i1.t2);
            const box2Current = this.computeIntervalBoundingBox(curve2, i2.t1, i2.t2);
            // Check if boxes intersect
            if (this.boundingBoxesIntersect(box1Current, box2Current)) {
                // Sample points along the curves to find intersection candidates
                const sampleCount = 20;
                const step1 = (i1.t2 - i1.t1) / sampleCount;
                const step2 = (i2.t2 - i2.t1) / sampleCount;
                for (let j = 0; j <= sampleCount; j++) {
                    const t1 = i1.t1 + j * step1;
                    const t2 = i2.t1 + j * step2;
                    const point1 = curve1.evaluate(t1);
                    const point2 = curve2.evaluate(t2);
                    // Check distance between points
                    const distance = point1.subtract(point2).length();
                    if (distance < this.tolerance) {
                        const pointKey = `${ point1.x.toFixed(6) },${ point1.y.toFixed(6) },${ point1.z.toFixed(6) }`;
                        if (!candidatePoints.has(pointKey)) {
                            candidatePoints.add(pointKey);
                            result.points.push(point1);
                            result.parameters.push({
                                curve1: t1,
                                curve2: t2
                            });
                        }
                    }
                }
            }
        }
        return result;
    }
}
class IntersectionResult {
    constructor() {
        // Initialize private storage with validation checks
        const validatePoints = points => {
            return Object.freeze(points.filter(point => {
                if (!(point instanceof Vector3D)) {
                    return false;
                }
                return Number.isFinite(point.x) && Number.isFinite(point.y) && Number.isFinite(point.z);
            }));
        };
        const validateParameters = params => {
            return Object.freeze(params.filter(param => {
                // For curve-curve intersections
                if (param.curve1 !== undefined && param.curve2 !== undefined) {
                    return Number.isFinite(param.curve1) && Number.isFinite(param.curve2);
                }
                // For surface-curve intersections
                if (param.curve !== undefined && param.surface !== undefined) {
                    return Number.isFinite(param.curve) && Number.isFinite(param.surface.u) && Number.isFinite(param.surface.v);
                }
                // For surface-surface intersections
                if (param.surface1 !== undefined && param.surface2 !== undefined) {
                    return Number.isFinite(param.surface1.u) && Number.isFinite(param.surface1.v) && Number.isFinite(param.surface2.u) && Number.isFinite(param.surface2.v);
                }
                return false;
            }));
        };
        const validateCurves = curves => {
            return Object.freeze(curves.filter(curve => {
                if (!Array.isArray(curve))
                    return false;
                return curve.every(segment => {
                    if (!segment.point || !segment.params)
                        return false;
                    if (!(segment.point instanceof Vector3D))
                        return false;
                    return Number.isFinite(segment.point.x) && Number.isFinite(segment.point.y) && Number.isFinite(segment.point.z);
                });
            }));
        };
        // Initialize properties with empty, immutable arrays
        Object.defineProperties(this, {
            points: {
                value: Object.freeze([]),
                writable: false,
                configurable: false
            },
            parameters: {
                value: Object.freeze([]),
                writable: false,
                configurable: false
            },
            curves: {
                value: Object.freeze([]),
                writable: false,
                configurable: false
            },
            validation: {
                value: {
                    validatePoints,
                    validateParameters,
                    validateCurves
                },
                writable: false,
                configurable: false
            },
            metadata: {
                value: Object.freeze({
                    createdAt: new Date(),
                    type: 'UNDEFINED',
                    precision: Number.EPSILON
                }),
                writable: false,
                configurable: false
            }
        });
        // Freeze the entire instance to prevent any modifications
        Object.freeze(this);
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
                        point: new Vector3D(segment.point.x === 0 ? 0 : segment.point.x, // Protect against -0
                        segment.point.y === 0 ? 0 : segment.point.y, segment.point.z === 0 ? 0 : segment.point.z),
                        params: Object.freeze({
                            surface1: { ...segment.params.surface1 },
                            surface2: { ...segment.params.surface2 }
                        })
                    };
                });
                // Validate all points in segment
                for (const segment of segmentPoints) {
                    if (!Number.isFinite(segment.point.x) || !Number.isFinite(segment.point.y) || !Number.isFinite(segment.point.z)) {
                        throw new Error('Non-finite values in intersection curve points');
                    }
                }
                // Return immutable curve segment
                return Object.freeze(segmentPoints);
            });
            // Return immutable array of intersection curves
            return Object.freeze(intersectionCurves);
        } catch (error) {
            throw new Error(`Failed to get intersection curves: ${ error.message }`);
        }
    }
    getParameters() {
        // If no parameters exist, return empty array
        if (!this.parameters || !Array.isArray(this.parameters)) {
            return Object.freeze([]);
        }
        try {
            // Deep copy and validate parameters
            const validatedParams = this.parameters.map(param => {
                // Validate parameter object exists
                if (!param || typeof param !== 'object') {
                    throw new Error('Invalid parameter object');
                }
                // Handle curve-curve intersection parameters
                if (param.curve1 !== undefined && param.curve2 !== undefined) {
                    if (!Number.isFinite(param.curve1) || !Number.isFinite(param.curve2)) {
                        throw new Error('Invalid curve-curve intersection parameters');
                    }
                    return Object.freeze({
                        curve1: param.curve1,
                        curve2: param.curve2
                    });
                }
                // Handle surface-curve intersection parameters
                if (param.curve !== undefined && param.surface !== undefined) {
                    if (!Number.isFinite(param.curve) || !param.surface || !Number.isFinite(param.surface.u) || !Number.isFinite(param.surface.v)) {
                        throw new Error('Invalid surface-curve intersection parameters');
                    }
                    return Object.freeze({
                        curve: param.curve,
                        surface: Object.freeze({
                            u: param.surface.u,
                            v: param.surface.v
                        })
                    });
                }
                // Handle surface-surface intersection parameters
                if (param.surface1 !== undefined && param.surface2 !== undefined) {
                    if (!param.surface1 || !param.surface2 || !Number.isFinite(param.surface1.u) || !Number.isFinite(param.surface1.v) || !Number.isFinite(param.surface2.u) || !Number.isFinite(param.surface2.v)) {
                        throw new Error('Invalid surface-surface intersection parameters');
                    }
                    return Object.freeze({
                        surface1: Object.freeze({
                            u: param.surface1.u,
                            v: param.surface1.v
                        }),
                        surface2: Object.freeze({
                            u: param.surface2.u,
                            v: param.surface2.v
                        })
                    });
                }
                throw new Error('Unrecognized parameter format');
            });
            // Return immutable array of validated parameters
            return Object.freeze(validatedParams);
        } catch (error) {
            throw new Error(`Failed to get intersection parameters: ${ error.message }`);
        }
    }
    classify() {
        try {
            // If no intersection data exists, return "NONE"
            if ((!this.points || this.points.length === 0) && (!this.curves || this.curves.length === 0)) {
                return {
                    type: 'NONE',
                    details: 'No intersection found'
                };
            }
            // If only discrete points exist
            if (this.points.length > 0 && (!this.curves || this.curves.length === 0)) {
                return {
                    type: 'POINT',
                    count: this.points.length,
                    details: `${ this.points.length } intersection point(s)`
                };
            }
            // If intersection curves exist
            if (this.curves && this.curves.length > 0) {
                // Check if curves form a closed loop
                const isClosedLoop = this.curves.some(curve => {
                    if (!Array.isArray(curve) || curve.length < 2)
                        return false;
                    const firstPoint = curve[0].point;
                    const lastPoint = curve[curve.length - 1].point;
                    return firstPoint.subtract(lastPoint).length() < Number.EPSILON;
                });
                if (isClosedLoop) {
                    return {
                        type: 'CLOSED_CURVE',
                        count: this.curves.length,
                        details: `${ this.curves.length } closed intersection curve(s)`
                    };
                } else {
                    return {
                        type: 'OPEN_CURVE',
                        count: this.curves.length,
                        details: `${ this.curves.length } open intersection curve(s)`
                    };
                }
            }
            // Check for tangential intersection
            if (this.parameters && this.parameters.length > 0) {
                const isTangential = this.parameters.some(param => {
                    // Check if derivatives are parallel at intersection point
                    if (param.surface1 && param.surface2) {
                        const derivatives1 = [
                            param.surface1.derivU,
                            param.surface1.derivV
                        ].filter(Boolean);
                        const derivatives2 = [
                            param.surface2.derivU,
                            param.surface2.derivV
                        ].filter(Boolean);
                        return derivatives1.some(d1 => derivatives2.some(d2 => Math.abs(Math.abs(d1.dot(d2)) - d1.length() * d2.length()) < Number.EPSILON));
                    }
                    return false;
                });
                if (isTangential) {
                    return {
                        type: 'TANGENTIAL',
                        count: this.points.length,
                        details: 'Tangential intersection'
                    };
                }
            }
            // If no specific classification can be determined
            return {
                type: 'UNKNOWN',
                details: 'Intersection type could not be determined'
            };
        } catch (error) {
            throw new Error(`Intersection classification failed: ${ error.message }`);
        }
    }
    // ... other methods ...
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
            throw new Error(`Failed to get intersection points: ${ error.message }`);
        }
    }
    getDistance(index = 0) {
        // Input validation
        if (!Number.isInteger(index) || index < 0) {
            throw new Error('Index must be a non-negative integer');
        }
        if (!this.points || index >= this.points.length) {
            throw new Error('Invalid intersection point index');
        }
        try {
            // Calculate distance between corresponding points
            if (this.points.length > index) {
                if (this.type === 'curve-curve') {
                    const curve1Point = this.points[index];
                    const curve2Param = this.parameters[index].curve2;
                    const curve2Point = this.curve2?.evaluate(curve2Param);
                    return curve1Point.subtract(curve2Point).length();
                } else if (this.type === 'surface-surface') {
                    const surface1Point = this.points[index];
                    const surface2Params = this.parameters[index].surface2;
                    const surface2Point = this.surface2?.evaluate(surface2Params.u, surface2Params.v);
                    return surface1Point.subtract(surface2Point).length();
                }
            }
            return 0;
        } catch (error) {
            throw new Error(`Failed to calculate intersection distance: ${ error.message }`);
        }
    }
    getAngle(index = 0) {
        // Input validation
        if (!Number.isInteger(index) || index < 0) {
            throw new Error('Index must be a non-negative integer');
        }
        if (!this.points || index >= this.points.length) {
            throw new Error('Invalid intersection point index');
        }
        try {
            // Calculate angle between intersecting entities
            if (this.points.length > index) {
                if (this.type === 'curve-curve') {
                    const t1 = this.parameters[index].curve1;
                    const t2 = this.parameters[index].curve2;
                    const tangent1 = this.curve1?.derivative(t1, 1);
                    const tangent2 = this.curve2?.derivative(t2, 1);
                    return Math.acos(Math.abs(tangent1.dot(tangent2)) / (tangent1.length() * tangent2.length()));
                } else if (this.type === 'surface-surface') {
                    const params1 = this.parameters[index].surface1;
                    const params2 = this.parameters[index].surface2;
                    const normal1 = this.surface1?.normal(params1.u, params1.v);
                    const normal2 = this.surface2?.normal(params2.u, params2.v);
                    return Math.acos(Math.abs(normal1.dot(normal2)) / (normal1.length() * normal2.length()));
                }
            }
            return 0;
        } catch (error) {
            throw new Error(`Failed to calculate intersection angle: ${ error.message }`);
        }
    }
    getMerged(tolerance = defaultTolerance) {
        // Merge nearby intersection points within tolerance
        if (!this.points || this.points.length === 0) {
            return this;
        }
        try {
            const mergedPoints = [];
            const mergedParams = [];
            const processed = new Set();
            for (let i = 0; i < this.points.length; i++) {
                if (processed.has(i))
                    continue;
                const point = this.points[i];
                const params = this.parameters[i];
                const group = {
                    points: [point],
                    params: [params]
                };
                // Find nearby points
                for (let j = i + 1; j < this.points.length; j++) {
                    if (processed.has(j))
                        continue;
                    if (point.subtract(this.points[j]).length() < tolerance) {
                        group.points.push(this.points[j]);
                        group.params.push(this.parameters[j]);
                        processed.add(j);
                    }
                }
                // Average the group
                if (group.points.length > 1) {
                    const avgPoint = group.points.reduce((acc, p) => acc.add(p), new Vector3D()).multiply(1 / group.points.length);
                    const avgParams = this.averageParameters(group.params);
                    mergedPoints.push(avgPoint);
                    mergedParams.push(avgParams);
                } else {
                    mergedPoints.push(point);
                    mergedParams.push(params);
                }
                processed.add(i);
            }
            // Create new result with merged data
            const mergedResult = new IntersectionResult();
            mergedResult.points = Object.freeze(mergedPoints);
            mergedResult.parameters = Object.freeze(mergedParams);
            mergedResult.type = this.type;
            return mergedResult;
        } catch (error) {
            throw new Error(`Failed to merge intersection points: ${ error.message }`);
        }
    }
    // Helper method for averaging parameters
    averageParameters(params) {
        if (!params || params.length === 0)
            return null;
        if (this.type === 'curve-curve') {
            return {
                curve1: params.reduce((acc, p) => acc + p.curve1, 0) / params.length,
                curve2: params.reduce((acc, p) => acc + p.curve2, 0) / params.length
            };
        } else if (this.type === 'surface-surface') {
            return {
                surface1: {
                    u: params.reduce((acc, p) => acc + p.surface1.u, 0) / params.length,
                    v: params.reduce((acc, p) => acc + p.surface1.v, 0) / params.length
                },
                surface2: {
                    u: params.reduce((acc, p) => acc + p.surface2.u, 0) / params.length,
                    v: params.reduce((acc, p) => acc + p.surface2.v, 0) / params.length
                }
            };
        }
        return null;
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
        // Check if knot vector exists
        if (!this.knots || this.knots.length === 0) {
            throw new Error('Cannot normalize empty knot vector');
        }
        try {
            // Get domain boundaries
            const firstKnot = this.knots[0];
            const lastKnot = this.knots[this.knots.length - 1];
            const domain = lastKnot - firstKnot;
            // Check for zero-length domain
            if (Math.abs(domain) < Number.EPSILON) {
                throw new Error('Cannot normalize knot vector with zero-length domain');
            }
            // Create new normalized knot vector
            const normalizedKnots = this.knots.map(knot => {
                // Normalize to [0,1] interval
                const normalizedValue = (knot - firstKnot) / domain;
                // Protect against numerical imprecision
                if (Math.abs(normalizedValue) < Number.EPSILON) {
                    return 0;
                }
                if (Math.abs(normalizedValue - 1) < Number.EPSILON) {
                    return 1;
                }
                return normalizedValue;
            });
            // Validate normalized knots maintain monotonic non-decreasing order
            for (let i = 1; i < normalizedKnots.length; i++) {
                if (normalizedKnots[i] < normalizedKnots[i - 1]) {
                    throw new Error('Normalization resulted in invalid knot sequence');
                }
            }
            // Return new KnotVector instance with normalized values
            return new KnotVector(normalizedKnots);
        } catch (error) {
            throw new Error(`Knot vector normalization failed: ${ error.message }`);
        }
    }
    insert(knot, multiplicity = 1) {
        // Input validation
        if (typeof knot !== 'number' || !Number.isFinite(knot)) {
            throw new Error('Knot must be a finite number');
        }
        if (!Number.isInteger(multiplicity) || multiplicity < 1) {
            throw new Error('Multiplicity must be a positive integer');
        }
        // Create a new knots array for the result
        const newKnots = [];
        let inserted = false;
        for (let i = 0; i < this.knots.length; i++) {
            // Insert the new knot at the correct position
            if (!inserted && (i === this.knots.length - 1 || this.knots[i] <= knot && this.knots[i + 1] >= knot)) {
                for (let j = 0; j < multiplicity; j++) {
                    newKnots.push(knot);
                }
                inserted = true;
            }
            // Mark that we have inserted the knot
            newKnots.push(this.knots[i]);
        }
        // Handle case where knot is inserted at the end
        if (!inserted) {
            for (let j = 0; j < multiplicity; j++) {
                newKnots.push(knot);
            }
        }
        // Create a new KnotVector with the updated knots
        return new KnotVector(newKnots);
    }
    remove(knot, multiplicity = 1) {
        // Input validation
        if (typeof knot !== 'number' || !Number.isFinite(knot)) {
            throw new Error('Knot must be a finite number');
        }
        if (!Number.isInteger(multiplicity) || multiplicity < 1) {
            throw new Error('Multiplicity must be a positive integer');
        }
        try {
            const knotsArray = [...this.knots];
            // Find indices of knots matching the input value within tolerance
            const matchingIndices = [];
            for (let i = 0; i < knotsArray.length; i++) {
                if (Math.abs(knotsArray[i] - knot) < Number.EPSILON) {
                    matchingIndices.push(i);
                }
            }
            // Check if we found enough occurrences to remove
            if (matchingIndices.length < multiplicity) {
                throw new Error('Cannot remove more knots than exist');
            }
            // Remove knots starting from the end to maintain index validity
            for (let i = 0; i < multiplicity; i++) {
                const indexToRemove = matchingIndices[matchingIndices.length - 1 - i];
                knotsArray.splice(indexToRemove, 1);
            }
            // Validate resulting knot vector
            this.validateKnotVector(knotsArray);
            // Create new KnotVector instance with updated knots
            return new KnotVector(knotsArray);
        } catch (error) {
            throw new Error(`Failed to remove knot: ${ error.message }`);
        }
    }
    refine(refinementPoints = [], minSpacing = 0.01) {
        // Input validation
        if (!Array.isArray(refinementPoints)) {
            throw new Error('Refinement points must be an array');
        }
        if (typeof minSpacing !== 'number' || minSpacing <= 0) {
            throw new Error('Minimum spacing must be a positive number');
        }
        try {
            // Filter and sort refinement points within domain
            const validPoints = refinementPoints.filter(point => {
                return typeof point === 'number' && Number.isFinite(point) && point >= this.knots[0] && point <= this.knots[this.knots.length - 1];
            }).sort((a, b) => a - b);
            // Initialize array for new knots
            const newKnots = [...this.knots];
            const insertedKnots = new Set();
            // Process each refinement point
            for (const point of validPoints) {
                // Find knot span containing the point
                let spanIndex = this.findKnotSpan(point);
                if (this.isTooCloseToExistingKnots(point, minSpacing)) {
                    continue;
                }
                // Skip if we've already inserted this point (within tolerance)
                const pointKey = point.toFixed(8);
                if (insertedKnots.has(pointKey)) {
                    continue;
                }
                // Insert new knot
                newKnots.splice(spanIndex + 1, 0, point);
                insertedKnots.add(pointKey);
            }
            // Add additional knots to maintain proper spacing
            this.insertIntermediateKnots(newKnots, minSpacing);
            // Ensure knots are properly ordered
            newKnots.sort((a, b) => a - b);
            // Create and return new knot vector
            return new KnotVector(newKnots);
        } catch (error) {
            throw new Error(`Knot vector refinement failed: ${ error.message }`);
        }
    }
    // ... other methods ...
    validate() {
        try {
            // Check if knots array exists and is not empty
            if (!this.knots || this.knots.length === 0) {
                throw new Error('Knot vector cannot be empty');
            }
            // Check minimum length (must be at least 2 knots for a valid knot vector)
            if (this.knots.length < 2) {
                throw new Error('Knot vector must contain at least 2 knots');
            }
            // Validate each knot value is a finite number
            for (let i = 0; i < this.knots.length; i++) {
                if (typeof this.knots[i] !== 'number' || !Number.isFinite(this.knots[i])) {
                    throw new Error(`Invalid knot value at index ${ i }: must be a finite number`);
                }
            }
            // Check for non-decreasing sequence
            for (let i = 1; i < this.knots.length; i++) {
                if (this.knots[i] < this.knots[i - 1]) {
                    throw new Error(`Invalid knot sequence at index ${ i }: knots must be non-decreasing`);
                }
            }
            // Check for valid domain
            const domain = this.domain;
            if (!Number.isFinite(domain.min) || !Number.isFinite(domain.max)) {
                throw new Error('Invalid knot vector domain');
            }
            if (domain.max <= domain.min) {
                throw new Error('Invalid domain: max must be greater than min');
            }
            // Check for valid multiplicity
            const multiplicities = this.multiplicities;
            if (!Array.isArray(multiplicities)) {
                throw new Error('Invalid multiplicities array');
            }
            // Calculate and validate sum of multiplicities
            const multiplicitySum = multiplicities.reduce((sum, mult) => sum + mult, 0);
            if (multiplicitySum !== this.knots.length) {
                throw new Error('Sum of multiplicities does not match knot vector length');
            }
            // Validate each multiplicity value
            for (let i = 0; i < multiplicities.length; i++) {
                if (!Number.isInteger(multiplicities[i]) || multiplicities[i] < 1) {
                    throw new Error(`Invalid multiplicity at index ${ i }: must be a positive integer`);
                }
            }
            // Check for no zero-length spans between unique knots
            let previousUniqueKnot = this.knots[0];
            for (let i = 1; i < this.knots.length; i++) {
                if (this.knots[i] > previousUniqueKnot) {
                    if (this.knots[i] - previousUniqueKnot < Number.EPSILON) {
                        throw new Error(`Invalid knot span at index ${ i }: zero-length span detected`);
                    }
                    previousUniqueKnot = this.knots[i];
                }
            }
            return {
                isValid: true,
                knots: this.knots.length,
                uniqueKnots: multiplicities.length,
                domain: { ...domain },
                multiplicities: [...multiplicities]
            };
        } catch (error) {
            return {
                isValid: false,
                error: error.message,
                knots: this.knots ? this.knots.length : 0,
                domain: this.domain ? { ...this.domain } : null
            };
        }
    }
    calculateMultiplicities() {
        // Input validation - ensure knots array exists
        if (!this.knots || !this.knots.length) {
            throw new Error('Cannot calculate multiplicities: knot vector is empty');
        }
        // Initialize array to store multiplicities
        const multiplicities = [];
        // Initialize tracking variables
        let currentKnot = this.knots[0];
        let currentMultiplicity = 1;
        // Iterate through knots starting from second knot
        for (let i = 1; i <= this.knots.length; i++) {
            if (i < this.knots.length && Math.abs(this.knots[i] - currentKnot) < Number.EPSILON) {
                // Current knot is same as previous (within tolerance)
                currentMultiplicity++;
            } else {
                // Different knot found or end reached
                multiplicities.push(currentMultiplicity);
                if (i < this.knots.length) {
                    // Start tracking new knot
                    currentKnot = this.knots[i];
                    currentMultiplicity = 1;
                }
            }
        }
        // Validate results
        const totalMultiplicity = multiplicities.reduce((sum, mult) => sum + mult, 0);
        if (totalMultiplicity !== this.knots.length) {
            throw new Error('Multiplicity calculation error: sum does not match knot vector length');
        }
        // Return immutable array of multiplicities
        return Object.freeze(multiplicities);
    }
    // ...existing methods...
    clone() {
        try {
            return new KnotVector([...this.knots]);
        } catch (error) {
            throw new Error(`Failed to clone knot vector: ${ error.message }`);
        }
    }
    getDomain() {
        try {
            return Object.freeze({ ...this.domain });
        } catch (error) {
            throw new Error(`Failed to get domain: ${ error.message }`);
        }
    }
    split(parameter, multiplicity = 1) {
        // Input validation
        if (typeof parameter !== 'number' || !Number.isFinite(parameter)) {
            throw new Error('Split parameter must be a finite number');
        }
        if (!Number.isInteger(multiplicity) || multiplicity < 1) {
            throw new Error('Multiplicity must be a positive integer');
        }
        // Check if parameter is within domain
        const min = this.knots[0];
        const max = this.knots[this.knots.length - 1];
        if (parameter <= min || parameter >= max) {
            throw new Error('Split parameter must be within knot vector domain (exclusive)');
        }
        try {
            // Find insertion index
            let insertionIndex = this.knots.findIndex(knot => knot > parameter);
            if (insertionIndex === -1) {
                throw new Error('Invalid split parameter');
            }
            // Calculate current multiplicity of parameter
            let currentMultiplicity = 0;
            for (const knot of this.knots) {
                if (Math.abs(knot - parameter) < Number.EPSILON) {
                    currentMultiplicity++;
                }
            }
            // Create new knot arrays for split result
            const leftKnots = [...this.knots.slice(0, insertionIndex)];
            const rightKnots = [...this.knots.slice(insertionIndex)];
            // Add parameter with specified multiplicity
            const knotsToAdd = multiplicity - currentMultiplicity;
            if (knotsToAdd > 0) {
                const insertKnots = new Array(knotsToAdd).fill(parameter);
                leftKnots.push(...insertKnots);
                rightKnots.unshift(...insertKnots);
            }
            // Create and return two new KnotVector instances
            return [
                new KnotVector(leftKnots),
                new KnotVector(rightKnots)
            ];
        } catch (error) {
            throw new Error(`Failed to split knot vector: ${ error.message }`);
        }
    }
    validateKnotVector(knots) {
        // Must have at least 2 knots
        if (knots.length < 2) {
            throw new Error('Knot vector must have at least 2 knots');
        }
        // Validate monotonic non-decreasing sequence
        for (let i = 1; i < knots.length; i++) {
            if (knots[i] < knots[i - 1]) {
                throw new Error('Knot vector must be monotonically non-decreasing');
            }
        }
        // Check for valid numerical values
        if (!knots.every(k => Number.isFinite(k))) {
            throw new Error('All knots must be finite numbers');
        }
        return true;
    }
    findKnotSpan(parameter) {
        // Find the knot span containing the parameter
        for (let i = 0; i < this.knots.length - 1; i++) {
            if (parameter >= this.knots[i] && parameter < this.knots[i + 1]) {
                return i;
            }
        }
        // If parameter equals last knot
        if (Math.abs(parameter - this.knots[this.knots.length - 1]) < Number.EPSILON) {
            return this.knots.length - 2;
        }
        return -1;
    }
    isTooCloseToExistingKnots(point, minSpacing) {
        return this.knots.some(knot => Math.abs(point - knot) < minSpacing);
    }
    insertIntermediateKnots(knots, minSpacing) {
        let i = 0;
        while (i < knots.length - 1) {
            const span = knots[i + 1] - knots[i];
            if (span > 2 * minSpacing) {
                // Calculate number of intermediate knots needed
                const numKnots = Math.floor(span / minSpacing) - 1;
                const spacing = span / (numKnots + 1);
                // Insert intermediate knots
                for (let j = 1; j <= numKnots; j++) {
                    const newKnot = knots[i] + j * spacing;
                    knots.splice(i + j, 0, newKnot);
                }
                i += numKnots + 1;
            } else {
                i++;
            }
        }
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
        try {
            const m = this.elements;
            const inv = new Float64Array(16);
            // Calculate adjoint matrix coefficients
            inv[0] = m[5] * m[10] * m[15] + m[6] * m[11] * m[13] + m[7] * m[9] * m[14] - m[7] * m[10] * m[13] - m[6] * m[9] * m[15] - m[5] * m[11] * m[14];
            inv[4] = m[4] * m[10] * m[15] + m[6] * m[11] * m[12] + m[7] * m[8] * m[14] - m[7] * m[10] * m[12] - m[6] * m[8] * m[15] - m[4] * m[11] * m[14];
            inv[8] = m[4] * m[9] * m[15] + m[5] * m[11] * m[12] + m[7] * m[8] * m[13] - m[7] * m[9] * m[12] - m[5] * m[8] * m[15] - m[4] * m[11] * m[13];
            inv[12] = m[4] * m[9] * m[14] + m[5] * m[10] * m[12] + m[6] * m[8] * m[13] - m[6] * m[9] * m[12] - m[5] * m[8] * m[14] - m[4] * m[10] * m[13];
            inv[1] = m[1] * m[10] * m[15] + m[2] * m[11] * m[13] + m[3] * m[9] * m[14] - m[3] * m[10] * m[13] - m[2] * m[9] * m[15] - m[1] * m[11] * m[14];
            inv[5] = m[0] * m[10] * m[15] + m[2] * m[11] * m[12] + m[3] * m[8] * m[14] - m[3] * m[10] * m[12] - m[2] * m[8] * m[15] - m[0] * m[11] * m[14];
            inv[9] = m[0] * m[9] * m[15] + m[1] * m[11] * m[12] + m[3] * m[8] * m[13] - m[3] * m[9] * m[12] - m[1] * m[8] * m[15] - m[0] * m[11] * m[13];
            inv[13] = m[0] * m[9] * m[14] + m[1] * m[10] * m[12] + m[2] * m[8] * m[13] - m[2] * m[9] * m[12] - m[1] * m[8] * m[14] - m[0] * m[10] * m[13];
            inv[2] = m[1] * m[6] * m[15] + m[2] * m[7] * m[13] + m[3] * m[5] * m[14] - m[3] * m[6] * m[13] - m[2] * m[5] * m[15] - m[1] * m[7] * m[14];
            inv[6] = m[0] * m[6] * m[15] + m[2] * m[7] * m[12] + m[3] * m[4] * m[14] - m[3] * m[6] * m[12] - m[2] * m[4] * m[15] - m[0] * m[7] * m[14];
            inv[10] = m[0] * m[5] * m[15] + m[1] * m[7] * m[12] + m[3] * m[4] * m[13] - m[3] * m[5] * m[12] - m[1] * m[4] * m[15] - m[0] * m[7] * m[13];
            inv[14] = m[0] * m[5] * m[14] + m[1] * m[6] * m[12] + m[2] * m[4] * m[13] - m[2] * m[5] * m[12] - m[1] * m[4] * m[14] - m[0] * m[6] * m[13];
            inv[3] = m[1] * m[6] * m[11] + m[2] * m[7] * m[9] + m[3] * m[5] * m[10] - m[3] * m[6] * m[9] - m[2] * m[5] * m[11] - m[1] * m[7] * m[10];
            inv[7] = m[0] * m[6] * m[11] + m[2] * m[7] * m[8] + m[3] * m[4] * m[10] - m[3] * m[6] * m[8] - m[2] * m[4] * m[11] - m[0] * m[7] * m[10];
            inv[11] = m[0] * m[5] * m[11] + m[1] * m[7] * m[8] + m[3] * m[4] * m[9] - m[3] * m[5] * m[8] - m[1] * m[4] * m[11] - m[0] * m[7] * m[9];
            inv[15] = m[0] * m[5] * m[10] + m[1] * m[6] * m[8] + m[2] * m[4] * m[9] - m[2] * m[5] * m[8] - m[1] * m[4] * m[10] - m[0] * m[6] * m[9];
            // Calculate determinant
            const det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
            // Check if matrix is invertible
            if (Math.abs(det) < Number.EPSILON) {
                throw new Error('Matrix is not invertible: determinant is zero');
            }
            // Calculate inverse by dividing by determinant
            const invDet = 1 / det;
            for (let i = 0; i < 16; i++) {
                inv[i] *= invDet;
                // Protect against -0
                if (Math.abs(inv[i]) < Number.EPSILON) {
                    inv[i] = 0;
                }
            }
            // Return new matrix with inverse elements
            return new Matrix4x4(inv);
        } catch (error) {
            throw new Error(`Matrix inversion failed: ${ error.message }`);
        }
    }
    transpose() {
        try {
            // Create new array for transposed elements
            const transposed = new Float64Array(16);
            // Transpose the matrix by swapping elements
            for (let i = 0; i < 4; i++) {
                for (let j = 0; j < 4; j++) {
                    // Convert from row-major to column-major order
                    transposed[j * 4 + i] = this.elements[i * 4 + j];
                }
            }
            // Create new Matrix4x4 instance with transposed elements
            return new Matrix4x4([
                transposed[0],
                transposed[1],
                transposed[2],
                transposed[3],
                transposed[4],
                transposed[5],
                transposed[6],
                transposed[7],
                transposed[8],
                transposed[9],
                transposed[10],
                transposed[11],
                transposed[12],
                transposed[13],
                transposed[14],
                transposed[15]
            ]);
        } catch (error) {
            throw new Error(`Matrix transposition failed: ${ error.message }`);
        }
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
            throw new Error(`Failed to create transformation matrix: ${ error.message }`);
        }
    }
    determinant() {
        try {
            const m = this.elements;
            // Calculate cofactors for first row expansion
            const coeff1 = m[0] * (m[5] * (m[10] * m[15] - m[11] * m[14]) - m[9] * (m[6] * m[15] - m[7] * m[14]) + m[13] * (m[6] * m[11] - m[7] * m[10]));
            const coeff2 = -m[4] * (m[1] * (m[10] * m[15] - m[11] * m[14]) - m[9] * (m[2] * m[15] - m[3] * m[14]) + m[13] * (m[2] * m[11] - m[3] * m[10]));
            const coeff3 = m[8] * (m[1] * (m[6] * m[15] - m[7] * m[14]) - m[5] * (m[2] * m[15] - m[3] * m[14]) + m[13] * (m[2] * m[7] - m[3] * m[6]));
            const coeff4 = -m[12] * (m[1] * (m[6] * m[11] - m[7] * m[10]) - m[5] * (m[2] * m[11] - m[3] * m[10]) + m[9] * (m[2] * m[7] - m[3] * m[6]));
            const det = coeff1 + coeff2 + coeff3 + coeff4;
            // Check for numerical validity
            if (!Number.isFinite(det)) {
                throw new Error('Determinant calculation resulted in non-finite value');
            }
            // Protect against -0
            return det === 0 ? 0 : det;
        } catch (error) {
            throw new Error(`Failed to calculate determinant: ${ error.message }`);
        }
    }
    // ...existing methods...
    clone() {
        try {
            return new Matrix4x4([...this.elements]);
        } catch (error) {
            throw new Error(`Failed to clone matrix: ${ error.message }`);
        }
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
                    throw new Error(`Invalid knot vector: ${ error.message }`);
                }
            } else {
                throw new Error('Knot vector must be an instance of KnotVector or an array');
            }
        }
        // Validate relationship between control points, degree, and knot vector
        const expectedKnotLength = controlPoints.length + degree + 1;
        if (knotVector.length !== expectedKnotLength) {
            throw new Error(`Knot vector length must be ${ expectedKnotLength } for ${ controlPoints.length } control points of degree ${ degree }`);
        }
        // Minimum number of control points must be degree + 1
        if (controlPoints.length < degree + 1) {
            throw new Error(`At least ${ degree + 1 } control points are required for degree ${ degree }`);
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
    stubMethod(){
        
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
            throw new Error(`Failed to create joined curve: ${ error.message }`);
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
            throw new Error(`Failed to create elevated curve: ${ error.message }`);
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
            throw new Error(`Failed to solve degree reduction system: ${ error.message }`);
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
            throw new Error(`Failed to create reduced curve: ${ error.message }`);
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
            throw new Error(`Failed to create reversed curve: ${ error.message }`);
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
    binomialCoefficient(n, k) {
        // Input validation
        if (!Number.isInteger(n) || !Number.isInteger(k)) {
            throw new Error('Binomial coefficient requires integer inputs');
        }
        if (k < 0 || k > n) {
            return 0;
        }
        // Optimize for edge cases
        if (k === 0 || k === n) {
            return 1;
        }
        if (k === 1 || k === n - 1) {
            return n;
        }
        // Optimize by using smaller k value due to symmetry
        k = Math.min(k, n - k);
        // Calculate using multiplicative formula
        let coefficient = 1;
        for (let i = 0; i < k; i++) {
            // Use intermediary calculation to prevent overflow
            coefficient *= n - i;
            coefficient /= i + 1;
            // Check for numerical overflow/underflow
            if (!Number.isFinite(coefficient)) {
                throw new Error('Binomial coefficient calculation resulted in non-finite value');
            }
        }
        // Handle potential floating point imprecision
        return Math.round(coefficient);
    }
    // ... other methods ...
    solveLU(A, b) {
        // Input validation
        if (!Array.isArray(A) || !Array.isArray(b)) {
            throw new Error('Matrix A and vector b must be arrays');
        }
        const n = A.length;
        // Validate matrix is square and matches b length
        if (!A.every(row => Array.isArray(row) && row.length === n)) {
            throw new Error('Matrix A must be square');
        }
        if (b.length !== n) {
            throw new Error('Vector b length must match matrix dimensions');
        }
        try {
            // Create augmented copies of arrays to prevent modifying input
            const LU = A.map(row => [...row]);
            const x = [...b];
            const y = new Array(n).fill(0);
            const pivot = new Array(n).fill(0);
            // Initialize pivot array
            for (let i = 0; i < n; i++) {
                pivot[i] = i;
            }
            // LU decomposition with partial pivoting
            for (let k = 0; k < n - 1; k++) {
                // Find pivot element
                let maxVal = Math.abs(LU[k][k]);
                let maxIndex = k;
                for (let i = k + 1; i < n; i++) {
                    const absVal = Math.abs(LU[i][k]);
                    if (absVal > maxVal) {
                        maxVal = absVal;
                        maxIndex = i;
                    }
                }
                // Check for singularity
                if (maxVal < Number.EPSILON) {
                    throw new Error('Matrix is singular or nearly singular');
                }
                // Swap rows if necessary
                if (maxIndex !== k) {
                    // Swap pivot indices
                    [pivot[k], pivot[maxIndex]] = [
                        pivot[maxIndex],
                        pivot[k]
                    ];
                    // Swap matrix rows
                    [LU[k], LU[maxIndex]] = [
                        LU[maxIndex],
                        LU[k]
                    ];
                }
                // Compute multipliers and eliminate k-th column
                for (let i = k + 1; i < n; i++) {
                    const multiplier = LU[i][k] / LU[k][k];
                    LU[i][k] = multiplier;
                    // Store multiplier in L part
                    // Update remaining elements in row i
                    for (let j = k + 1; j < n; j++) {
                        LU[i][j] -= multiplier * LU[k][j];
                    }
                }
            }
            // Forward substitution to solve Ly = b
            for (let i = 0; i < n; i++) {
                y[i] = x[pivot[i]];
                for (let j = 0; j < i; j++) {
                    y[i] -= LU[i][j] * y[j];
                }
            }
            // Back substitution to solve Ux = y
            for (let i = n - 1; i >= 0; i--) {
                for (let j = i + 1; j < n; j++) {
                    y[i] -= LU[i][j] * y[j];
                }
                // Check for division by zero
                if (Math.abs(LU[i][i]) < Number.EPSILON) {
                    throw new Error('Division by zero in back substitution');
                }
                y[i] /= LU[i][i];
            }
            // Ensure all values are finite
            if (!y.every(val => Number.isFinite(val))) {
                throw new Error('Solution contains non-finite values');
            }
            return y;
        } catch (error) {
            throw new Error(`LU decomposition failed: ${ error.message }`);
        }
    }
    // ...existing methods...
    length() {
        try {
            // Calculate curve length using adaptive Gaussian quadrature
            const gaussPoints = [
                {
                    t: -0.861136311594053,
                    w: 0.347854845137454
                },
                {
                    t: -0.339981043584856,
                    w: 0.652145154862546
                },
                {
                    t: 0.339981043584856,
                    w: 0.652145154862546
                },
                {
                    t: 0.861136311594053,
                    w: 0.347854845137454
                }
            ];
            const numSegments = Math.ceil((this.domain.max - this.domain.min) * 10);
            const dt = (this.domain.max - this.domain.min) / numSegments;
            let length = 0;
            for (let i = 0; i < numSegments; i++) {
                const t0 = this.domain.min + i * dt;
                const t1 = t0 + dt;
                const tmid = (t0 + t1) / 2;
                const hdt = dt / 2;
                let segmentLength = 0;
                for (const gp of gaussPoints) {
                    const t = tmid + hdt * gp.t;
                    const derivative = this.derivative(t, 1);
                    segmentLength += gp.w * derivative.length();
                }
                length += segmentLength * hdt;
            }
            return length === 0 ? 0 : length;
        } catch (error) {
            throw new Error(`Failed to calculate curve length: ${ error.message }`);
        }
    }
    isPointOnCurve(point, tolerance = defaultTolerance) {
        try {
            if (!(point instanceof Vector3D)) {
                throw new Error('Point must be a Vector3D instance');
            }
            const closest = this.closestPoint(point, tolerance);
            return closest.distance <= tolerance;
        } catch (error) {
            throw new Error(`Failed to check if point is on curve: ${ error.message }`);
        }
    }
    enforceContinuity(continuity, direction) {
        // Input validation
        if (!Number.isInteger(continuity) || continuity < 0 || continuity > 2) {
            throw new Error('Continuity must be 0 (C0), 1 (C1), or 2 (C2)');
        }
        if (direction !== 'start' && direction !== 'end') {
            throw new Error('Direction must be either "start" or "end"');
        }
        try {
            // Get affected control points based on continuity
            const numPoints = this.controlPoints.length;
            const pointsToModify = continuity + 1;
            if (pointsToModify >= numPoints) {
                throw new Error('Not enough control points to enforce requested continuity');
            }
            // Get reference points and derivatives
            const t = direction === 'start' ? this.domain.min : this.domain.max;
            const point = this.evaluate(t);
            const derivatives = [];
            for (let i = 1; i <= continuity; i++) {
                derivatives.push(this.derivative(t, i));
            }
            // Modify control points to enforce continuity
            if (direction === 'start') {
                // Enforce continuity at start
                for (let i = 0; i <= continuity; i++) {
                    if (i === 0) {
                        // C0 continuity: position
                        const [x, y, z] = point;
                        const w = this.controlPoints[0].weight();
                        this.controlPoints[0] = new ControlPoint(x, y, z, w);
                    } else {
                        // C1 and C2 continuity: derivatives
                        const derivative = derivatives[i - 1];
                        const scale = 1 / factorial(i);
                        const baseWeight = this.controlPoints[i].weight();
                        // Calculate new control point position based on derivative
                        const newPosition = point.add(derivative.multiply(scale));
                        this.controlPoints[i] = new ControlPoint(newPosition.x, newPosition.y, newPosition.z, baseWeight);
                    }
                }
            } else {
                // Enforce continuity at end
                for (let i = 0; i <= continuity; i++) {
                    const index = numPoints - 1 - i;
                    if (i === 0) {
                        // C0 continuity: position
                        const [x, y, z] = point;
                        const w = this.controlPoints[index].weight();
                        this.controlPoints[index] = new ControlPoint(x, y, z, w);
                    } else {
                        // C1 and C2 continuity: derivatives
                        const derivative = derivatives[i - 1];
                        const scale = Math.pow(-1, i) / factorial(i);
                        const baseWeight = this.controlPoints[index].weight();
                        // Calculate new control point position based on derivative
                        const newPosition = point.add(derivative.multiply(scale));
                        this.controlPoints[index] = new ControlPoint(newPosition.x, newPosition.y, newPosition.z, baseWeight);
                    }
                }
            }
            // Return the modified curve
            return new NurbsCurve(this.controlPoints, this.knotVector, this.degree);
        } catch (error) {
            throw new Error(`Failed to enforce continuity: ${ error.message }`);
        }
    }
    execute(parameter) {
        // Input validation
        if (typeof parameter !== 'number' || !Number.isFinite(parameter)) {
            throw new Error('Parameter must be a finite number');
        }
        // Check if curve is properly initialized
        if (!this.controlPoints || !this.knotVector || !Number.isInteger(this.degree)) {
            throw new Error('Curve is not properly initialized');
        }
        try {
            // Determine operation based on parameter range
            if (parameter <= this.domain.min) {
                // Return start point
                return this.evaluate(this.domain.min);
            }
            if (parameter >= this.domain.max) {
                // Return end point
                return this.evaluate(this.domain.max);
            }
            // For parameters within domain
            // Calculate point on curve
            const point = this.evaluate(parameter);
            // Calculate first derivative (tangent)
            const tangent = this.derivative(parameter, 1);
            // Calculate second derivative (curvature)
            const secondDeriv = this.derivative(parameter, 2);
            // Calculate curvature vector
            let curvatureVector = null;
            const tangentLength = tangent.length();
            if (tangentLength > Number.EPSILON) {
                const crossProduct = tangent.cross(secondDeriv);
                const curvature = crossProduct.length() / Math.pow(tangentLength, 3);
                if (curvature > Number.EPSILON) {
                    curvatureVector = crossProduct.normalize().multiply(curvature);
                }
            }
            // Return execution result object
            return Object.freeze({
                parameter: parameter,
                point: point,
                tangent: tangent.normalize(),
                curvature: curvatureVector,
                success: true
            });
        } catch (error) {
            throw new Error(`Curve execution failed: ${ error.message }`);
        }
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
            throw new Error(`Failed to calculate surface normal: ${ error.message }`);
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
            throw new Error(`Failed to create split surfaces: ${ error.message }`);
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
            throw new Error(`Failed to calculate surface curvature: ${ error.message }`);
        }
    }
    // ...existing methods...
    clone() {
        try {
            const clonedSurface = new NurbsSurface();
            clonedSurface.controlPoints = this.controlPoints.map(row => row.map(cp => cp.clone()));
            clonedSurface.knotVectorU = this.knotVectorU.clone();
            clonedSurface.knotVectorV = this.knotVectorV.clone();
            clonedSurface.degreeU = this.degreeU;
            clonedSurface.degreeV = this.degreeV;
            return clonedSurface;
        } catch (error) {
            throw new Error(`Failed to clone surface: ${ error.message }`);
        }
    }
    isPointOnSurface(point, tolerance = defaultTolerance) {
        try {
            if (!(point instanceof Vector3D)) {
                throw new Error('Point must be a Vector3D instance');
            }
            // Sample parameters to find approximate closest point
            const numSamples = 10;
            const uStep = (this.knotVectorU.domain.max - this.knotVectorU.domain.min) / numSamples;
            const vStep = (this.knotVectorV.domain.max - this.knotVectorV.domain.min) / numSamples;
            let minDistance = Infinity;
            for (let i = 0; i <= numSamples; i++) {
                for (let j = 0; j <= numSamples; j++) {
                    const u = this.knotVectorU.domain.min + i * uStep;
                    const v = this.knotVectorV.domain.min + j * vStep;
                    const surfacePoint = this.evaluate(u, v);
                    const distance = point.subtract(surfacePoint).length();
                    minDistance = Math.min(minDistance, distance);
                }
            }
            return minDistance <= tolerance;
        } catch (error) {
            throw new Error(`Failed to check if point is on surface: ${ error.message }`);
        }
    }
    join(otherSurface, direction = 'u', continuity = 1) {
        // Input validation
        if (!(otherSurface instanceof NurbsSurface)) {
            throw new Error('Join operation requires another NurbsSurface instance');
        }
        if (direction !== 'u' && direction !== 'v') {
            throw new Error('Direction must be either "u" or "v"');
        }
        if (!Number.isInteger(continuity) || continuity < 0 || continuity > 2) {
            throw new Error('Continuity must be 0, 1, or 2');
        }
        try {
            // Check if surfaces are compatible for joining
            if (direction === 'u') {
                // For U-direction join, V dimensions must match
                if (this.controlPoints[0].length !== otherSurface.controlPoints[0].length || this.degreeV !== otherSurface.degreeV) {
                    throw new Error('Surfaces are not compatible for U-direction joining');
                }
                // Get end points of this surface and start points of other surface
                const thisEnd = this.evaluate(this.knotVectorU.domain.max, this.knotVectorV.domain.min);
                const otherStart = otherSurface.evaluate(otherSurface.knotVectorU.domain.min, otherSurface.knotVectorV.domain.min);
                // Check if surfaces share boundary (within tolerance)
                if (thisEnd.subtract(otherStart).length() > defaultTolerance) {
                    throw new Error('Surfaces must share a common boundary for joining');
                }
            } else {
                // For V-direction join, U dimensions must match
                if (this.controlPoints.length !== otherSurface.controlPoints.length || this.degreeU !== otherSurface.degreeU) {
                    throw new Error('Surfaces are not compatible for V-direction joining');
                }
                // Get end points of this surface and start points of other surface
                const thisEnd = this.evaluate(this.knotVectorU.domain.min, this.knotVectorV.domain.max);
                const otherStart = otherSurface.evaluate(otherSurface.knotVectorU.domain.min, otherSurface.knotVectorV.domain.min);
                // Check if surfaces share boundary (within tolerance)
                if (thisEnd.subtract(otherStart).length() > defaultTolerance) {
                    throw new Error('Surfaces must share a common boundary for joining');
                }
            }
            // Create new control points array
            let newControlPoints = [];
            let newKnotVectorU, newKnotVectorV;
            if (direction === 'u') {
                // Merge control points in U direction
                newControlPoints = [...this.controlPoints];
                // Remove last row if ensuring continuity
                if (continuity > 0)
                    newControlPoints.pop();
                newControlPoints.push(...otherSurface.controlPoints);
                // Merge U knot vectors
                const firstMax = this.knotVectorU.domain.max;
                const secondMin = otherSurface.knotVectorU.domain.min;
                const secondScale = otherSurface.knotVectorU.domain.max - secondMin;
                const adjustedSecondKnots = otherSurface.knotVectorU.knots.map(knot => {
                    const normalized = (knot - secondMin) / secondScale;
                    return firstMax + normalized * secondScale;
                });
                newKnotVectorU = new KnotVector([
                    ...this.knotVectorU.knots,
                    ...adjustedSecondKnots.slice(this.degreeU + 1)
                ]);
                newKnotVectorV = this.knotVectorV;
            } else {
                // Merge control points in V direction
                newControlPoints = this.controlPoints.map((row, i) => {
                    const newRow = [...row];
                    // Remove last column if ensuring continuity
                    if (continuity > 0)
                        newRow.pop();
                    return [
                        ...newRow,
                        ...otherSurface.controlPoints[i]
                    ];
                });
                // Keep U knot vector, merge V knot vectors
                const firstMax = this.knotVectorV.domain.max;
                const secondMin = otherSurface.knotVectorV.domain.min;
                const secondScale = otherSurface.knotVectorV.domain.max - secondMin;
                const adjustedSecondKnots = otherSurface.knotVectorV.knots.map(knot => {
                    const normalized = (knot - secondMin) / secondScale;
                    return firstMax + normalized * secondScale;
                });
                newKnotVectorU = this.knotVectorU;
                newKnotVectorV = new KnotVector([
                    ...this.knotVectorV.knots,
                    ...adjustedSecondKnots.slice(this.degreeV + 1)
                ]);
            }
            // Create joined surface
            const joinedSurface = new NurbsSurface();
            joinedSurface.controlPoints = newControlPoints;
            joinedSurface.knotVectorU = newKnotVectorU;
            joinedSurface.knotVectorV = newKnotVectorV;
            joinedSurface.degreeU = this.degreeU;
            joinedSurface.degreeV = this.degreeV;
            // Ensure continuity at join
            if (continuity > 0) {
                joinedSurface.enforceContinuity(continuity, direction);
            }
            return joinedSurface;
        } catch (error) {
            throw new Error(`Surface joining failed: ${ error.message }`);
        }
    }
    reverse(direction = 'both') {
        // Input validation
        if (![
                'u',
                'v',
                'both'
            ].includes(direction)) {
            throw new Error('Direction must be "u", "v", or "both"');
        }
        try {
            // Create new surface with same degrees
            const reversedSurface = new NurbsSurface();
            reversedSurface.degreeU = this.degreeU;
            reversedSurface.degreeV = this.degreeV;
            // Helper function to reverse knot vector
            const reverseKnotVector = knotVector => {
                const knots = knotVector.knots;
                const a = knots[0];
                const b = knots[knots.length - 1];
                const reversedKnots = new Float64Array(knots.length);
                for (let i = 0; i < knots.length; i++) {
                    reversedKnots[i] = a + b - knots[knots.length - 1 - i];
                }
                return new KnotVector(reversedKnots);
            };
            // Reverse U direction if needed
            if (direction === 'u' || direction === 'both') {
                reversedSurface.knotVectorU = reverseKnotVector(this.knotVectorU);
            } else {
                reversedSurface.knotVectorU = this.knotVectorU.clone();
            }
            // Reverse V direction if needed
            if (direction === 'v' || direction === 'both') {
                reversedSurface.knotVectorV = reverseKnotVector(this.knotVectorV);
            } else {
                reversedSurface.knotVectorV = this.knotVectorV.clone();
            }
            // Create reversed control point grid
            const numU = this.controlPoints.length;
            const numV = this.controlPoints[0].length;
            const reversedPoints = Array(numU).fill(null).map(() => Array(numV));
            // Copy and reverse control points based on direction
            for (let i = 0; i < numU; i++) {
                for (let j = 0; j < numV; j++) {
                    let sourceI = i;
                    let sourceJ = j;
                    if (direction === 'u' || direction === 'both') {
                        sourceI = numU - 1 - i;
                    }
                    if (direction === 'v' || direction === 'both') {
                        sourceJ = numV - 1 - j;
                    }
                    reversedPoints[i][j] = this.controlPoints[sourceI][sourceJ].clone();
                }
            }
            reversedSurface.controlPoints = reversedPoints;
            // Validate the reversed surface
            const validationResult = this.validateSurfaceData(reversedSurface.controlPoints, reversedSurface.knotVectorU, reversedSurface.knotVectorV, reversedSurface.degreeU, reversedSurface.degreeV);
            if (!validationResult) {
                throw new Error('Reversed surface validation failed');
            }
            return reversedSurface;
        } catch (error) {
            throw new Error(`Surface reversal failed: ${ error.message }`);
        }
    }
    // Helper method to validate surface data
    validateSurfaceData(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV) {
        if (!Array.isArray(controlPoints) || controlPoints.length === 0) {
            return false;
        }
        const numU = controlPoints.length;
        const numV = controlPoints[0].length;
        // Check control point grid is rectangular
        if (!controlPoints.every(row => Array.isArray(row) && row.length === numV)) {
            return false;
        }
        // Check control points are valid
        if (!controlPoints.every(row => row.every(point => point instanceof ControlPoint))) {
            return false;
        }
        // Check knot vector compatibility
        if (knotVectorU.length !== numU + degreeU + 1) {
            return false;
        }
        if (knotVectorV.length !== numV + degreeV + 1) {
            return false;
        }
        // Check degrees are valid
        if (!Number.isInteger(degreeU) || degreeU < 1 || !Number.isInteger(degreeV) || degreeV < 1) {
            return false;
        }
        return true;
    }
    enforceContinuity(continuity, direction) {
        // Input validation
        if (!Number.isInteger(continuity) || continuity < 0 || continuity > 2) {
            throw new Error('Continuity must be 0 (C0), 1 (C1), or 2 (C2)');
        }
        if (direction !== 'u' && direction !== 'v') {
            throw new Error('Direction must be either "u" or "v"');
        }
        try {
            // Determine number of control points needed based on continuity
            const numPointsToModify = continuity + 1;
            if (direction === 'u') {
                // Get reference points at start/end of U direction
                const startU = this.knotVectorU.knots[this.degreeU];
                const endU = this.knotVectorU.knots[this.knotVectorU.length - this.degreeU - 1];
                // Calculate derivatives at U boundaries for each V parameter
                for (let j = 0; j < this.controlPoints[0].length; j++) {
                    // Start boundary (u = startU)
                    const startDerivatives = [];
                    for (let k = 0; k <= continuity; k++) {
                        startDerivatives.push(this.derivative(startU, this.knotVectorV.knots[j], k, 0));
                    }
                    // End boundary (u = endU)
                    const endDerivatives = [];
                    for (let k = 0; k <= continuity; k++) {
                        endDerivatives.push(this.derivative(endU, this.knotVectorV.knots[j], k, 0));
                    }
                    // Modify control points to enforce continuity
                    for (let k = 0; k < numPointsToModify; k++) {
                        // Start boundary control points
                        const scale = 1 / factorial(k);
                        const startNewPoint = startDerivatives[k].multiply(scale);
                        const baseWeight = this.controlPoints[k][j].weight();
                        this.controlPoints[k][j] = new ControlPoint(startNewPoint.x, startNewPoint.y, startNewPoint.z, baseWeight);
                        // End boundary control points
                        const endScale = Math.pow(-1, k) / factorial(k);
                        const endNewPoint = endDerivatives[k].multiply(endScale);
                        const endWeight = this.controlPoints[this.controlPoints.length - 1 - k][j].weight();
                        this.controlPoints[this.controlPoints.length - 1 - k][j] = new ControlPoint(endNewPoint.x, endNewPoint.y, endNewPoint.z, endWeight);
                    }
                }
            } else {
                // Get reference points at start/end of V direction
                const startV = this.knotVectorV.knots[this.degreeV];
                const endV = this.knotVectorV.knots[this.knotVectorV.length - this.degreeV - 1];
                // Calculate derivatives at V boundaries for each U parameter
                for (let i = 0; i < this.controlPoints.length; i++) {
                    // Start boundary (v = startV)
                    const startDerivatives = [];
                    for (let k = 0; k <= continuity; k++) {
                        startDerivatives.push(this.derivative(this.knotVectorU.knots[i], startV, 0, k));
                    }
                    // End boundary (v = endV)
                    const endDerivatives = [];
                    for (let k = 0; k <= continuity; k++) {
                        endDerivatives.push(this.derivative(this.knotVectorU.knots[i], endV, 0, k));
                    }
                    // Modify control points to enforce continuity
                    for (let k = 0; k < numPointsToModify; k++) {
                        // Start boundary control points
                        const scale = 1 / factorial(k);
                        const startNewPoint = startDerivatives[k].multiply(scale);
                        const baseWeight = this.controlPoints[i][k].weight();
                        this.controlPoints[i][k] = new ControlPoint(startNewPoint.x, startNewPoint.y, startNewPoint.z, baseWeight);
                        // End boundary control points
                        const endScale = Math.pow(-1, k) / factorial(k);
                        const endNewPoint = endDerivatives[k].multiply(endScale);
                        const endWeight = this.controlPoints[i][this.controlPoints[i].length - 1 - k].weight();
                        this.controlPoints[i][this.controlPoints[i].length - 1 - k] = new ControlPoint(endNewPoint.x, endNewPoint.y, endNewPoint.z, endWeight);
                    }
                }
            }
            // Return modified surface
            return this;
        } catch (error) {
            throw new Error(`Failed to enforce continuity: ${ error.message }`);
        }
    }
    // ... other properties and methods ...
    tessellate(tolerance = defaultTolerance) {
        // Input validation
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        const points = [];
        const parameters = [];
        const triangles = [];
        const uRange = {
            min: this.knotVectorU.knots[this.degreeU],
            max: this.knotVectorU.knots[this.knotVectorU.length - this.degreeU - 1]
        };
        const vRange = {
            min: this.knotVectorV.knots[this.degreeV],
            max: this.knotVectorV.knots[this.knotVectorV.length - this.degreeV - 1]
        };
        const subdividePatch = (u1, u2, v1, v2, depth = 0) => {
            if (depth > 10)
                return;
            // Prevent infinite recursion
            const p00 = this.evaluate(u1, v1);
            const p10 = this.evaluate(u2, v1);
            const p01 = this.evaluate(u1, v2);
            const p11 = this.evaluate(u2, v2);
            const uMid = (u1 + u2) / 2;
            const vMid = (v1 + v2) / 2;
            const pMidMid = this.evaluate(uMid, vMid);
            const interpolated = p00.add(p10).add(p01).add(p11).multiply(0.25);
            const deviation = pMidMid.subtract(interpolated).length();
            if (deviation > tolerance) {
                subdividePatch(u1, uMid, v1, vMid, depth + 1);
                subdividePatch(uMid, u2, v1, vMid, depth + 1);
                subdividePatch(u1, uMid, vMid, v2, depth + 1);
                subdividePatch(uMid, u2, vMid, v2, depth + 1);
            } else {
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
        subdividePatch(uRange.min, uRange.max, vRange.min, vRange.max);
        return {
            points: Object.freeze(points),
            parameters: Object.freeze(parameters),
            triangles: Object.freeze(triangles)
        };
    }
}
class Shell {
    constructor() {
        this.faces = [];
        this.edges = [];
        this.vertices = [];
    }
    addFace(face) {
        if (!(face instanceof TopologyFace)) {
            throw new Error('Parameter must be an instance of TopologyFace');
        }
        // Ensure the face does not already exist in the shell
        if (this.faces.includes(face)) {
            throw new Error('Face already exists in the shell');
        }
        // Add the face to the shell's collection
        this.faces.push(face);
        // Update edges and vertices collections based on the face's edges
        face.edges.forEach(edge => {
            if (!this.edges.includes(edge)) {
                this.edges.push(edge);
            }
            edge.faces.push(face);
        });
        // Ensure the edge recognizes the face
        face.bounds.forEach(loop => {
            loop.edges.forEach(edge => {
                if (!this.edges.includes(edge)) {
                    this.edges.push(edge);
                }
            });
        });
        // Update vertices collection based on edges
        this.vertices = Array.from(new Set([
            ...this.vertices,
            ...face.edges.flatMap(edge => [
                edge.startVertex,
                edge.endVertex
            ])
        ]));
    }
    removeFace(face) {
        // Input validation
        if (!(face instanceof TopologyFace)) {
            throw new Error('Parameter must be an instance of TopologyFace');
        }
        const index = this.faces.indexOf(face);
        // Check if the face exists in the shell
        if (index === -1) {
            throw new Error('The specified face does not exist in the shell');
        }
        // Remove the face from the shell's collection
        this.faces.splice(index, 1);
        // Remove edges associated with the face
        face.edges.forEach(edge => {
            const edgeIndex = this.edges.indexOf(edge);
            if (edgeIndex !== -1) {
                this.edges.splice(edgeIndex, 1);
            }
            // Remove face reference from the edge
            const faceIndex = edge.faces.indexOf(face);
            if (faceIndex !== -1) {
                edge.faces.splice(faceIndex, 1);
            }
        });
        // Update the vertices collection
        this.vertices = Array.from(new Set([
            ...this.vertices,
            ...face.edges.flatMap(edge => [
                edge.startVertex,
                edge.endVertex
            ])
        ]));
    }
    heal() {
        // Input validation
        if (this.faces.length === 0) {
            throw new Error('Cannot heal: no faces in the shell');
        }
        // Initialize a set to track the removed edges
        const edgesToRemove = new Set();
        // Iterate over faces to find gaps between edges
        for (const face of this.faces) {
            const faceEdges = face.edges;
            for (let i = 0; i < faceEdges.length; i++) {
                const currentEdge = faceEdges[i];
                const nextEdge = faceEdges[(i + 1) % faceEdges.length];
                // Check for gaps between current and next edge
                const currentEnd = currentEdge.endVertex.position;
                const nextStart = nextEdge.startVertex.position;
                if (currentEnd.subtract(nextStart).length() > defaultTolerance) {
                    // Attempt to heal gap by repositioning vertices
                    const midPoint = new Vector3D((currentEnd.x + nextStart.x) / 2, (currentEnd.y + nextStart.y) / 2, (currentEnd.z + nextStart.z) / 2);
                    // Update the next edge's start vertex
                    nextEdge.startVertex.setPosition(midPoint);
                }
            }
        }
        // Clean up any duplicate edges and vertices
        this.edges = [...new Set(this.edges)];
        this.vertices = [...new Set(this.vertices)];
        // Validate final state after healing
        this.validate();
    }
    // ... other properties and methods ...
    validate() {
        if (this.faces.length === 0) {
            throw new Error('Shell must contain at least one face');
        }
        // Check that each face has associated edges
        for (const face of this.faces) {
            if (!face.edges || face.edges.length === 0) {
                throw new Error('Each face must have at least one associated edge');
            }
        }
        // Check that vertices are connected through edges
        const vertexSet = new Set(this.vertices);
        if (vertexSet.size === 0) {
            throw new Error('Shell must contain at least one vertex');
        }
        for (const edge of this.edges) {
            if (!vertexSet.has(edge.startVertex) || !vertexSet.has(edge.endVertex)) {
                throw new Error('All edges must connect valid vertices in the shell');
            }
        }
        // Check for duplicate edges
        const uniqueEdges = new Set(this.edges);
        if (uniqueEdges.size !== this.edges.length) {
            throw new Error('Duplicate edges detected in the shell');
        }
        // Check for manifoldness: Each edge should belong to exactly two faces
        const edgeFaceMap = new Map();
        for (const edge of this.edges) {
            const faceCount = edge.faces.length;
            if (faceCount > 2) {
                throw new Error('Non-manifold edge detected: an edge can belong to at most two faces');
            }
            edgeFaceMap.set(edge, faceCount);
        }
        return true;
    }
}
class Solid {
    constructor() {
        this.shells = [];
    }
    validate() {
        if (this.shells.length === 0) {
            throw new Error('Solid must contain at least one shell');
        }
    }
    heal() {
        if (this.shells.length === 0) {
            throw new Error('Cannot heal: no shells in the solid');
        }
        let modified = false;
        // Iterate through each shell to heal it
        for (const shell of this.shells) {
            try {
                shell.heal();
                modified = true;
            } // Mark as modified if any shell was healed
            catch (error) {
                console.warn(`Failed to heal shell: ${ error.message }`);
            }
        }
        // If any modifications were made, validate the solid
        if (modified) {
            this.validate();
        }
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
            throw new Error(`Tessellation failed: ${ error.message }`);
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
            throw new Error(`Surface tessellation failed: ${ error.message }`);
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
            throw new Error(`Tessellation of trimmed surface failed: ${ error.message }`);
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
            throw new Error(`Failed to generate knot lines: ${ error.message }`);
        }
    }
}
class TopologyEdge {
    constructor() {
        this.curve = null;
        this.startVertex = null;
        this.endVertex = null;
        this.faces = [];
    }
    reverse() {
        try {
            // Validate edge state
            if (!this.curve || !this.startVertex || !this.endVertex) {
                throw new Error('Edge is not properly initialized');
            }
            // Reverse the NURBS curve
            this.curve = this.curve.reverse();
            // Swap vertices
            const tempVertex = this.startVertex;
            this.startVertex = this.endVertex;
            this.endVertex = tempVertex;
            // Update vertex edge references
            this.startVertex.edges.forEach((edge, index, edges) => {
                if (edge === this) {
                    edges[index] = this;
                }
            });
            this.endVertex.edges.forEach((edge, index, edges) => {
                if (edge === this) {
                    edges[index] = this;
                }
            });
            // Update face orientation references if needed
            this.faces.forEach(face => {
                if (face.bounds) {
                    face.bounds.forEach(loop => {
                        if (loop.edges) {
                            const edgeIndex = loop.edges.indexOf(this);
                            if (edgeIndex !== -1) {
                                loop.edges[edgeIndex] = this;
                            }
                        }
                    });
                }
            });
            // Validate new state
            if (!this.validateReversedState()) {
                throw new Error('Invalid edge state after reversal');
            }
            return this;
        } catch (error) {
            throw new Error(`Edge reversal failed: ${ error.message }`);
        }
    }
    split(parameter) {
        // Input validation
        if (typeof parameter !== 'number' || !Number.isFinite(parameter)) {
            throw new Error('Split parameter must be a finite number');
        }
        // Validate edge state
        if (!this.curve || !this.startVertex || !this.endVertex) {
            throw new Error('Edge is not properly initialized');
        }
        // Check if parameter is within curve domain
        if (parameter <= this.curve.domain.min || parameter >= this.curve.domain.max) {
            throw new Error('Split parameter must be within curve domain (exclusive of endpoints)');
        }
        try {
            // Split the underlying NURBS curve
            const [curve1, curve2] = this.curve.split(parameter);
            // Create new vertex at split point
            const splitPoint = this.curve.evaluate(parameter);
            const newVertex = new TopologyVertex({
                position: splitPoint,
                tolerance: Math.min(this.startVertex.tolerance, this.endVertex.tolerance)
            });
            // Create new edges
            const edge1 = new TopologyEdge();
            const edge2 = new TopologyEdge();
            // Set curves
            edge1.curve = curve1;
            edge2.curve = curve2;
            // Set vertices for first edge
            edge1.startVertex = this.startVertex;
            edge1.endVertex = newVertex;
            // Set vertices for second edge
            edge2.startVertex = newVertex;
            edge2.endVertex = this.endVertex;
            // Update vertex-edge references
            newVertex.edges = [
                edge1,
                edge2
            ];
            // Update start vertex edge reference
            const startEdgeIndex = this.startVertex.edges.indexOf(this);
            if (startEdgeIndex !== -1) {
                this.startVertex.edges[startEdgeIndex] = edge1;
            }
            // Update end vertex edge reference
            const endEdgeIndex = this.endVertex.edges.indexOf(this);
            if (endEdgeIndex !== -1) {
                this.endVertex.edges[endEdgeIndex] = edge2;
            }
            // Handle face connections
            edge1.faces = [...this.faces];
            edge2.faces = [...this.faces];
            // Update face edge references
            this.faces.forEach(face => {
                if (face.bounds) {
                    face.bounds.forEach(loop => {
                        if (loop.edges) {
                            const edgeIndex = loop.edges.indexOf(this);
                            if (edgeIndex !== -1) {
                                // Replace current edge with two new edges in the loop
                                loop.edges.splice(edgeIndex, 1, edge1, edge2);
                            }
                        }
                    });
                }
            });
            // Validate new topology
            if (!this.validateSplitResult(edge1, edge2, newVertex)) {
                throw new Error('Invalid topology after split');
            }
            // Clear references from original edge
            this.curve = null;
            this.startVertex = null;
            this.endVertex = null;
            this.faces = [];
            // Return new elements
            return {
                edge1,
                edge2,
                vertex: newVertex
            };
        } catch (error) {
            throw new Error(`Edge split failed: ${ error.message }`);
        }
    }
    merge(otherEdge) {
        // Input validation
        if (!(otherEdge instanceof TopologyEdge)) {
            throw new Error('Parameter must be a TopologyEdge instance');
        }
        if (this.startVertex !== otherEdge.startVertex && this.endVertex !== otherEdge.startVertex && this.startVertex !== otherEdge.endVertex && this.endVertex !== otherEdge.endVertex) {
            throw new Error('Edges must share a common vertex to merge');
        }
        // Determine merging direction
        let newCurve = this.curve.join(otherEdge.curve);
        let newStartVertex = this.startVertex;
        let newEndVertex = this.endVertex;
        // If they don't share a start vertex, adjust the new start and end vertices
        if (this.endVertex === otherEdge.startVertex) {
            newEndVertex = otherEdge.endVertex;
        } else if (this.startVertex === otherEdge.endVertex) {
            newStartVertex = otherEdge.startVertex;
        } else if (this.endVertex === otherEdge.endVertex) {
            newEndVertex = this.startVertex;
            newStartVertex = otherEdge.startVertex;
        }
        // Create a new merged edge
        const mergedEdge = new TopologyEdge();
        mergedEdge.curve = newCurve;
        mergedEdge.startVertex = newStartVertex;
        mergedEdge.endVertex = newEndVertex;
        // Merge faces
        mergedEdge.faces = Array.from(new Set([
            ...this.faces,
            ...otherEdge.faces
        ]));
        // Update the edge references in start and end vertices
        newStartVertex.edges.push(mergedEdge);
        newEndVertex.edges.push(mergedEdge);
        // Invalidate original edges
        this.curve = null;
        this.startVertex = null;
        this.endVertex = null;
        this.faces = [];
        otherEdge.curve = null;
        otherEdge.startVertex = null;
        otherEdge.endVertex = null;
        otherEdge.faces = [];
        return mergedEdge;
    }
    getFaces() {
        // Return a frozen copy of the faces array to prevent modification
        return Object.freeze([...this.faces]);
    }
    validateReversedState() {
        // Verify curve endpoints match vertex positions
        const startParam = this.curve.domain.min;
        const endParam = this.curve.domain.max;
        const startPoint = this.curve.evaluate(startParam);
        const endPoint = this.curve.evaluate(endParam);
        const tolerance = 1e-7;
        // Use appropriate tolerance
        // Check start vertex position matches curve start
        if (startPoint.subtract(this.startVertex.position).length() > tolerance) {
            return false;
        }
        // Check end vertex position matches curve end
        if (endPoint.subtract(this.endVertex.position).length() > tolerance) {
            return false;
        }
        // Verify vertices reference this edge
        if (!this.startVertex.edges.includes(this) || !this.endVertex.edges.includes(this)) {
            return false;
        }
        // Verify face references are maintained
        if (!this.faces.every(face => face.edges.includes(this))) {
            return false;
        }
        return true;
    }
    // Helper method to validate split result
    validateSplitResult(edge1, edge2, vertex) {
        const tolerance = vertex.tolerance;
        // Verify curve continuity at split point
        const endParam1 = edge1.curve.domain.max;
        const startParam2 = edge2.curve.domain.min;
        const endPoint1 = edge1.curve.evaluate(endParam1);
        const startPoint2 = edge2.curve.evaluate(startParam2);
        // Check position continuity
        if (endPoint1.subtract(startPoint2).length() > tolerance) {
            return false;
        }
        // Check tangent continuity
        const tangent1 = edge1.curve.derivative(endParam1, 1).normalize();
        const tangent2 = edge2.curve.derivative(startParam2, 1).normalize();
        if (Math.abs(Math.abs(tangent1.dot(tangent2)) - 1) > tolerance) {
            return false;
        }
        // Verify vertex connections
        if (!edge1.startVertex.edges.includes(edge1) || !edge1.endVertex.edges.includes(edge1) || !edge2.startVertex.edges.includes(edge2) || !edge2.endVertex.edges.includes(edge2)) {
            return false;
        }
        // Verify face connections
        if (!edge1.faces.every(face => face.edges.includes(edge1)) || !edge2.faces.every(face => face.edges.includes(edge2))) {
            return false;
        }
        return true;
    }
    mergeEdges(edge1, edge2) {
        // Input validation
        if (!(edge1 instanceof TopologyEdge) || !(edge2 instanceof TopologyEdge)) {
            throw new Error('Both parameters must be TopologyEdge instances');
        }
        // Check if edges share a common vertex
        const commonVertex = this.findCommonVertex(edge1, edge2);
        if (!commonVertex) {
            throw new Error('Edges must share a common vertex to merge');
        }
        try {
            // Create a new edge for the merged result
            const mergedEdge = new TopologyEdge();
            // Determine merge direction and set vertices
            if (edge1.endVertex === commonVertex && edge2.startVertex === commonVertex) {
                mergedEdge.startVertex = edge1.startVertex;
                mergedEdge.endVertex = edge2.endVertex;
                mergedEdge.curve = edge1.curve.join(edge2.curve);
            } else if (edge1.startVertex === commonVertex && edge2.endVertex === commonVertex) {
                mergedEdge.startVertex = edge2.startVertex;
                mergedEdge.endVertex = edge1.endVertex;
                mergedEdge.curve = edge2.curve.join(edge1.curve);
            } else if (edge1.startVertex === commonVertex && edge2.startVertex === commonVertex) {
                edge2.reverse();
                mergedEdge.startVertex = edge2.endVertex;
                mergedEdge.endVertex = edge1.endVertex;
                mergedEdge.curve = edge2.curve.join(edge1.curve);
            } else if (edge1.endVertex === commonVertex && edge2.endVertex === commonVertex) {
                edge2.reverse();
                mergedEdge.startVertex = edge1.startVertex;
                mergedEdge.endVertex = edge2.startVertex;
                mergedEdge.curve = edge1.curve.join(edge2.curve);
            } else {
                throw new Error('Invalid edge configuration for merging');
            }
            // Update vertex references
            mergedEdge.startVertex.edges = mergedEdge.startVertex.edges.filter(e => e !== edge1 && e !== edge2);
            mergedEdge.startVertex.edges.push(mergedEdge);
            mergedEdge.endVertex.edges = mergedEdge.endVertex.edges.filter(e => e !== edge1 && e !== edge2);
            mergedEdge.endVertex.edges.push(mergedEdge);
            // Merge face references
            mergedEdge.faces = Array.from(new Set([
                ...edge1.faces,
                ...edge2.faces
            ]));
            // Update face references
            mergedEdge.faces.forEach(face => {
                const index1 = face.edges.indexOf(edge1);
                const index2 = face.edges.indexOf(edge2);
                if (index1 !== -1)
                    face.edges[index1] = mergedEdge;
                if (index2 !== -1)
                    face.edges.splice(index2, 1);
            });
            // Invalidate original edges
            edge1.curve = null;
            edge1.startVertex = null;
            edge1.endVertex = null;
            edge1.faces = [];
            edge2.curve = null;
            edge2.startVertex = null;
            edge2.endVertex = null;
            edge2.faces = [];
            return mergedEdge;
        } catch (error) {
            throw new Error(`Edge merge failed: ${ error.message }`);
        }
    }
    findSharedVertex(edge1, edge2) {
        if (edge1.startVertex === edge2.startVertex || edge1.startVertex === edge2.endVertex) {
            return edge1.startVertex;
        }
        if (edge1.endVertex === edge2.startVertex || edge1.endVertex === edge2.endVertex) {
            return edge1.endVertex;
        }
        return null;
    }
    validateMergedEdge(edge) {
        // Verify edge has valid curve
        if (!edge.curve || !(edge.curve instanceof NurbsCurve)) {
            throw new Error('Invalid curve in merged edge');
        }
        // Verify vertex connections
        if (!edge.startVertex || !edge.endVertex) {
            throw new Error('Invalid vertex connections in merged edge');
        }
        // Verify edge is in vertex edge lists
        if (!edge.startVertex.edges.includes(edge) || !edge.endVertex.edges.includes(edge)) {
            throw new Error('Edge not properly referenced in vertex lists');
        }
        // Verify face references
        if (!edge.faces.every(face => face.edges.includes(edge))) {
            throw new Error('Invalid face references in merged edge');
        }
    }
    findCommonVertex(edge1, edge2) {
        if (edge1.startVertex === edge2.startVertex || edge1.startVertex === edge2.endVertex) {
            return edge1.startVertex;
        }
        if (edge1.endVertex === edge2.startVertex || edge1.endVertex === edge2.endVertex) {
            return edge1.endVertex;
        }
        return null;
    }
}
class TopologyFace {
    constructor() {
        this.surface = null;
        this.edges = [];
        this.bounds = [];
    }
    addBound(loop) {
        if (!(loop instanceof TopologyLoop)) {
            throw new Error('Parameter must be an instance of TopologyLoop');
        }
        this.bounds.push(loop);
        this.validateBounds();
    }
    removeBound(loop) {
        if (!(loop instanceof TopologyLoop)) {
            throw new Error('Parameter must be an instance of TopologyLoop');
        }
        const index = this.bounds.indexOf(loop);
        if (index === -1) {
            throw new Error('The specified loop is not associated with this face');
        }
        this.bounds.splice(index, 1);
        this.validateBounds();
    }
    split(criteria) {
        // Input validation
        if (typeof criteria !== 'function') {
            throw new Error('Criteria must be a function');
        }
        // Array to hold new faces created from the split
        const newFaces = [];
        // Iterate over edges to determine potential split points
        this.edges.forEach(edge => {
            const splitPoints = criteria(edge);
            // Use criteria function to get split points
            if (splitPoints) {
                // Create new faces based on split points
                const newFace = new TopologyFace();
                // Assign edges based on the split point
                this.edges.forEach(e => {
                    if (splitPoints.includes(e)) {
                        newFace.addBound(e);
                    } else {
                        this.addBound(e);
                    }
                });
                newFaces.push(newFace);
            }
        });
        // Validate and return new faces
        if (newFaces.length === 0) {
            throw new Error('No valid splits were performed');
        }
        return newFaces;
    }
    merge(face) {
        // Input validation
        if (!(face instanceof TopologyFace)) {
            throw new Error('Parameter must be a TopologyFace instance');
        }
        // Ensure faces have the same surface
        if (this.surface !== face.surface) {
            throw new Error('Faces must belong to the same surface to merge');
        }
        // Combine edges
        const combinedEdges = new Set([
            ...this.edges,
            ...face.edges
        ]);
        this.edges = Array.from(combinedEdges);
        // Combine bounds
        const combinedBounds = new Set([
            ...this.bounds,
            ...face.bounds
        ]);
        this.bounds = Array.from(combinedBounds);
        // Validate the merged face
        this.validateBounds();
    }
    validateBounds() {
        const uniqueLoops = new Set(this.bounds);
        if (uniqueLoops.size !== this.bounds.length) {
            throw new Error('Duplicate bounds found in face');
        }
    }
}
class TopologyLoop {
    constructor() {
        this.edges = [];
        this.face = null;
    }
    reverse() {
        // Validate the loop has edges
        if (this.edges.length === 0) {
            throw new Error('Cannot reverse an empty loop');
        }
        // Reverse the order of edges in the loop
        this.edges.reverse();
        // Additionally, reverse the orientation of each edge
        this.edges.forEach(edge => edge.reverse());
    }
    split(criteria) {
        // Input validation
        if (typeof criteria !== 'function') {
            throw new Error('Criteria must be a function');
        }
        const newLoops = [];
        let currentLoopEdges = [];
        for (const edge of this.edges) {
            if (criteria(edge)) {
                // If criteria is met, create a new loop for the split
                if (currentLoopEdges.length > 0) {
                    newLoops.push(new TopologyLoop());
                    newLoops[newLoops.length - 1].edges = currentLoopEdges;
                    currentLoopEdges = [];
                }
            }
            currentLoopEdges.push(edge);
        }
        // If there are any remaining edges, create a final loop
        if (currentLoopEdges.length > 0) {
            const finalLoop = new TopologyLoop();
            finalLoop.edges = currentLoopEdges;
            newLoops.push(finalLoop);
        }
        // Update this loop's edges to empty and validate the new loops
        this.edges = [];
        newLoops.forEach(loop => loop.validateEdges());
        return newLoops;
    }
    merge(loop) {
        // Input validation
        if (!(loop instanceof TopologyLoop)) {
            throw new Error('Parameter must be an instance of TopologyLoop');
        }
        // Ensure loops have the same face
        if (this.face !== loop.face) {
            throw new Error('Cannot merge loops that belong to different faces');
        }
        try {
            // Find connecting edges between the loops
            const connectingEdges = this.edges.filter(edge1 => loop.edges.some(edge2 => edge1.endVertex === edge2.startVertex || edge1.endVertex === edge2.endVertex || edge1.startVertex === edge2.startVertex || edge1.startVertex === edge2.endVertex));
            if (connectingEdges.length === 0) {
                throw new Error('Loops must share at least one vertex to merge');
            }
            // Create ordered set of edges for the merged loop
            const mergedEdges = [...this.edges];
            const remainingEdges = [...loop.edges];
            // Remove duplicate edges
            remainingEdges.forEach(edge => {
                if (!mergedEdges.some(existing => existing.startVertex === edge.startVertex && existing.endVertex === edge.endVertex)) {
                    mergedEdges.push(edge);
                }
            });
            // Update edge connectivity
            for (let i = 0; i < mergedEdges.length; i++) {
                const currentEdge = mergedEdges[i];
                const nextEdge = mergedEdges[(i + 1) % mergedEdges.length];
                // Ensure edges are connected
                if (currentEdge.endVertex !== nextEdge.startVertex) {
                    // If edges aren't connected, try to reverse the next edge
                    if (currentEdge.endVertex === nextEdge.endVertex) {
                        nextEdge.reverse();
                    } else {
                        throw new Error('Cannot create continuous loop from merged edges');
                    }
                }
            }
            // Update this loop's edges
            this.edges = mergedEdges;
            // Clear the other loop's edges
            loop.edges = [];
            loop.face = null;
            // Validate the merged loop
            this.validateEdges();
            return this;
        } catch (error) {
            throw new Error(`Loop merge failed: ${ error.message }`);
        }
    }
    validateEdges() {
        // Ensure that all edges are valid and belong to this face
        if (!this.edges.length) {
            throw new Error('Loop must contain at least one edge');
        }
        const uniqueEdges = new Set(this.edges);
        if (uniqueEdges.size !== this.edges.length) {
            throw new Error('Duplicate edges found in loop');
        }
    }
}
class TopologyVertex {
    constructor(options = {}) {
        // Extract options with defaults
        const {position = null, edges = [], tolerance = 0.000001, id = crypto.randomUUID()} = options;
        // Validate position if provided
        if (position !== null && !(position instanceof Vector3D)) {
            throw new Error('Position must be a Vector3D instance');
        }
        // Validate edges array
        if (!Array.isArray(edges)) {
            throw new Error('Edges must be an array');
        }
        if (!edges.every(edge => edge instanceof TopologyEdge)) {
            throw new Error('All edges must be TopologyEdge instances');
        }
        // Validate tolerance
        if (typeof tolerance !== 'number' || tolerance <= 0) {
            throw new Error('Tolerance must be a positive number');
        }
        // Initialize private properties with validation
        Object.defineProperties(this, {
            '_position': {
                value: position ? position.clone() : null,
                writable: true,
                enumerable: false
            },
            '_edges': {
                value: [...edges],
                writable: true,
                enumerable: false
            },
            '_tolerance': {
                value: tolerance,
                writable: false,
                enumerable: false
            },
            '_id': {
                value: id,
                writable: false,
                enumerable: false
            },
            '_timestamp': {
                value: Date.now(),
                writable: false,
                enumerable: false
            }
        });
        // Validate initial state
        this.validateState();
    }
    merge(otherVertex) {
        // Input validation
        if (!(otherVertex instanceof TopologyVertex)) {
            throw new Error('Parameter must be a TopologyVertex instance');
        }
        if (!this.position || !otherVertex.position) {
            throw new Error('Both vertices must have valid positions');
        }
        // Check if vertices are close enough to merge
        const distance = this.position.subtract(otherVertex.position).length();
        if (distance > this.tolerance) {
            throw new Error('Vertices are too far apart to merge');
        }
        try {
            // Store original states for rollback if needed
            const originalPosition = this.position.clone();
            const originalEdges = [...this.edges];
            const otherOriginalEdges = [...otherVertex.edges];
            // Calculate new position based on edge curvature and tangent continuity
            const newPosition = this.calculateOptimalMergePosition(otherVertex);
            this.position = newPosition;
            // Track processed edges to avoid duplicates
            const processedEdges = new Set();
            // First, handle edges that will be merged
            const edgePairs = this.findMergeableEdgePairs(otherVertex);
            for (const [edge1, edge2] of edgePairs) {
                this.mergeEdges(edge1, edge2);
                processedEdges.add(edge1);
                processedEdges.add(edge2);
            }
            // Update remaining edges from both vertices
            for (const edge of [
                    ...this.edges,
                    ...otherVertex.edges
                ]) {
                if (!processedEdges.has(edge)) {
                    if (edge.startVertex === otherVertex) {
                        edge.startVertex = this;
                    }
                    if (edge.endVertex === otherVertex) {
                        edge.endVertex = this;
                    }
                    // Adjust edge geometry if needed
                    this.adjustEdgeGeometry(edge);
                    if (!this.edges.includes(edge)) {
                        this.edges.push(edge);
                    }
                }
            }
            // Sort edges by parametric value for consistency
            this.sortEdgesByParameter();
            // Validate merged state
            try {
                this.validateMergedTopology();
            } catch (error) {
                // Rollback if validation fails
                this.position = originalPosition;
                this.edges = originalEdges;
                otherVertex.edges = otherOriginalEdges;
                throw error;
            }
            // Clear other vertex's connections
            otherVertex.edges = [];
            otherVertex.position = null;
            return this;
        } catch (error) {
            throw new Error(`Vertex merge failed: ${ error.message }`);
        }
    }
    split() {
        // Validate current state
        if (!this.position || this._edges.length < 2) {
            throw new Error('Vertex must have a valid position and at least two edges to split');
        }
        try {
            // Create new vertex with same position and tolerance
            const newVertex = new TopologyVertex({
                position: this.position.clone(),
                tolerance: this.tolerance
            });
            // Analyze edge tangent directions to determine optimal split
            const edgeGroups = this.groupEdgesByDirection();
            // Validate we have at least two distinct groups
            if (edgeGroups.length < 2) {
                throw new Error('Cannot split vertex: edges do not form distinct groups');
            }
            // Keep first group with this vertex, move second group to new vertex
            const [group1, group2] = edgeGroups;
            // Update edge connections for second group to new vertex
            for (const edge of group2) {
                if (edge.startVertex === this) {
                    edge.startVertex = newVertex;
                }
                if (edge.endVertex === this) {
                    edge.endVertex = newVertex;
                }
            }
            // Update edge lists
            this._edges = group1;
            newVertex._edges = group2;
            // Validate topology after split
            this.validateState();
            newVertex.validateState();
            // Return both vertices
            return {
                originalVertex: this,
                newVertex: newVertex
            };
        } catch (error) {
            throw new Error(`Vertex split failed: ${ error.message }`);
        }
    }
    getEdges(options = {}) {
        // Input validation
        if (typeof options !== 'object') {
            throw new Error('Options parameter must be an object');
        }
        try {
            // Create copy of edges array for manipulation
            let edges = [...this._edges];
            // Filter edges if filter function provided
            if (typeof options.filter === 'function') {
                edges = edges.filter(options.filter);
            }
            // Handle different sorting options
            if (options.sortBy) {
                switch (options.sortBy) {
                case 'parameter':
                    // Sort by curve parameter
                    edges.sort((a, b) => {
                        const paramA = a.startVertex === this ? a.curve.domain.min : a.curve.domain.max;
                        const paramB = b.startVertex === this ? b.curve.domain.min : b.curve.domain.max;
                        return paramA - paramB;
                    });
                    break;
                case 'angle':
                    // Sort by angle from reference direction
                    if (!options.referenceDirection || !(options.referenceDirection instanceof Vector3D)) {
                        throw new Error('Reference direction must be provided for angle sorting');
                    }
                    const reference = options.referenceDirection.normalize();
                    edges.sort((a, b) => {
                        const angleA = this.getEdgeAngle(a, reference);
                        const angleB = this.getEdgeAngle(b, reference);
                        return angleA - angleB;
                    });
                    break;
                case 'length':
                    // Sort by edge curve length
                    edges.sort((a, b) => {
                        return a.curve.length() - b.curve.length();
                    });
                    break;
                default:
                    throw new Error(`Invalid sort option: ${ options.sortBy }`);
                }
                // Apply reverse sorting if specified
                if (options.reverse) {
                    edges.reverse();
                }
            }
            // Additional metadata if requested
            if (options.includeMetadata) {
                return edges.map(edge => ({
                    edge: edge,
                    isStart: edge.startVertex === this,
                    tangent: this.getEdgeTangent(edge),
                    parameter: edge.startVertex === this ? edge.curve.domain.min : edge.curve.domain.max,
                    faces: edge.faces
                }));
            }
            // Return the processed edges array
            return Object.freeze(edges);
        } catch (error) {
            throw new Error(`Failed to get edges: ${ error.message }`);
        }
    }
    setPosition(newPosition) {
        // Input validation
        if (!(newPosition instanceof Vector3D)) {
            throw new Error('New position must be a Vector3D instance');
        }
        if (!Number.isFinite(newPosition.x) || !Number.isFinite(newPosition.y) || !Number.isFinite(newPosition.z)) {
            throw new Error('Position coordinates must be finite numbers');
        }
        try {
            // Store original position for rollback if needed
            const originalPosition = this._position ? this._position.clone() : null;
            // Update position
            this._position = newPosition.clone();
            // Check if movement is within tolerance for connected edges
            for (const edge of this._edges) {
                // Get curve parameter for this vertex
                const parameter = edge.startVertex === this ? edge.curve.domain.min : edge.curve.domain.max;
                // Get point on curve at parameter
                const curvePoint = edge.curve.evaluate(parameter);
                // Check distance between new position and curve endpoint
                const distance = curvePoint.subtract(newPosition).length();
                if (distance > this._tolerance) {
                    // If movement exceeds tolerance, attempt to adjust connected edges
                    if (!this.adjustConnectedEdges(edge, newPosition)) {
                        // Rollback position if adjustment fails
                        this._position = originalPosition;
                        throw new Error('Cannot maintain edge connectivity with new position');
                    }
                }
            }
            // Validate new state after position change
            this.validateState();
            // Update any dependent geometry
            this.updateDependentGeometry();
            return true;
        } catch (error) {
            throw new Error(`Failed to set position: ${ error.message }`);
        }
    }
    validateMergedTopology() {
        // Check for invalid edge connections
        for (const edge of this.edges) {
            if (edge.startVertex !== this && edge.endVertex !== this) {
                throw new Error('Invalid edge connection after merge');
            }
        }
        // Check for invalid edge intersections
        for (let i = 0; i < this.edges.length; i++) {
            for (let j = i + 1; j < this.edges.length; j++) {
                const curve1 = this.edges[i].curve;
                const curve2 = this.edges[j].curve;
                // Skip if curves are the same
                if (curve1 === curve2)
                    continue;
                // Check for tangency at vertex
                const param1 = this.edges[i].startVertex === this ? curve1.domain.min : curve1.domain.max;
                const param2 = this.edges[j].startVertex === this ? curve2.domain.min : curve2.domain.max;
                const tangent1 = curve1.derivative(param1, 1).normalize();
                const tangent2 = curve2.derivative(param2, 1).normalize();
                const angleDeviation = Math.abs(Math.abs(tangent1.dot(tangent2)) - 1);
                if (angleDeviation < this.tolerance) {
                    throw new Error('Invalid tangent edge configuration after merge');
                }
            }
        }
        // Validate position is finite
        if (!this.position || !Number.isFinite(this.position.x) || !Number.isFinite(this.position.y) || !Number.isFinite(this.position.z)) {
            throw new Error('Invalid vertex position after merge');
        }
    }
    // Public getter for position
    get position() {
        return this._position ? this._position.clone() : null;
    }
    // Public getter for edges
    get edges() {
        return [...this._edges];
    }
    // Public getter for tolerance
    get tolerance() {
        return this._tolerance;
    }
    // Public getter for id
    get id() {
        return this._id;
    }
    // Helper method to validate vertex state
    validateState() {
        // Check position-edge consistency
        if (this._position && this._edges.length > 0) {
            // Verify all connected edges actually connect to this vertex
            const invalidEdges = this._edges.filter(edge => edge.startVertex !== this && edge.endVertex !== this);
            if (invalidEdges.length > 0) {
                throw new Error('Invalid edge connections detected');
            }
            // Verify edge endpoints are within tolerance of vertex position
            for (const edge of this._edges) {
                const edgePoint = edge.startVertex === this ? edge.curve.evaluate(edge.curve.domain.min) : edge.curve.evaluate(edge.curve.domain.max);
                const distance = edgePoint.subtract(this._position).length();
                if (distance > this._tolerance) {
                    throw new Error('Edge endpoint deviation exceeds tolerance');
                }
            }
        }
        // Check for duplicate edges
        const edgeSet = new Set(this._edges);
        if (edgeSet.size !== this._edges.length) {
            throw new Error('Duplicate edges detected');
        }
        // Additional geometric validation
        if (this._position) {
            if (!Number.isFinite(this._position.x) || !Number.isFinite(this._position.y) || !Number.isFinite(this._position.z)) {
                throw new Error('Vertex position contains non-finite coordinates');
            }
        }
        return true;
    }
    // Helper method to group edges by direction
    groupEdgesByDirection() {
        // Initialize groups array
        const groups = [];
        const processedEdges = new Set();
        // Process each edge
        for (const edge of this._edges) {
            if (processedEdges.has(edge))
                continue;
            // Get tangent vector at vertex
            const tangent = this.getEdgeTangent(edge);
            const currentGroup = [edge];
            processedEdges.add(edge);
            // Find other edges with similar direction
            for (const otherEdge of this._edges) {
                if (processedEdges.has(otherEdge))
                    continue;
                const otherTangent = this.getEdgeTangent(otherEdge);
                // Check if tangents are parallel or anti-parallel
                const dotProduct = Math.abs(tangent.dot(otherTangent));
                if (Math.abs(dotProduct - 1) < this.tolerance) {
                    currentGroup.push(otherEdge);
                    processedEdges.add(otherEdge);
                }
            }
            groups.push(currentGroup);
        }
        // Sort groups by size (descending) for consistency
        groups.sort((a, b) => b.length - a.length);
        return groups;
    }
    // Helper method to get edge tangent at vertex
    getEdgeTangent(edge) {
        try {
            const isStart = edge.startVertex === this;
            const parameter = isStart ? edge.curve.domain.min : edge.curve.domain.max;
            const tangent = edge.curve.derivative(parameter, 1);
            // Normalize and handle direction based on vertex position
            const normalized = tangent.normalize();
            return isStart ? normalized : normalized.multiply(-1);
        } catch (error) {
            throw new Error(`Failed to calculate edge tangent: ${ error.message }`);
        }
    }
    // Helper method to calculate angle between edge and reference direction
    getEdgeAngle(edge, referenceDirection) {
        const tangent = this.getEdgeTangent(edge);
        const angle = Math.acos(Math.min(1, Math.max(-1, tangent.dot(referenceDirection))));
        // Determine sign of angle using cross product
        const cross = referenceDirection.cross(tangent);
        const sign = cross.z >= 0 ? 1 : -1;
        return sign * angle;
    }
    // Helper method to adjust connected edges when position changes
    adjustConnectedEdges(edge, newPosition) {
        try {
            const isStart = edge.startVertex === this;
            const parameter = isStart ? edge.curve.domain.min : edge.curve.domain.max;
            // Get current curve type
            if (edge.curve instanceof NurbsCurve) {
                // For NURBS curves, adjust the nearest control point
                const controlPoints = [...edge.curve.controlPoints];
                const index = isStart ? 0 : controlPoints.length - 1;
                // Create new control point at new position
                const weight = controlPoints[index].weight();
                controlPoints[index] = new ControlPoint(newPosition.x, newPosition.y, newPosition.z, weight);
                // Create new curve with adjusted control point
                const newCurve = new NurbsCurve(controlPoints, edge.curve.knotVector, edge.curve.degree);
                // Verify new curve maintains continuity with adjacent geometry
                if (this.validateCurveContinuity(newCurve, edge)) {
                    edge.curve = newCurve;
                    return true;
                }
            }
            // If adjustment fails or curve type is unsupported
            return false;
        } catch (error) {
            return false;
        }
    }
    // Helper method to validate curve continuity
    validateCurveContinuity(newCurve, edge) {
        try {
            // Check geometric continuity with adjacent edges
            for (const adjacentEdge of this._edges) {
                if (adjacentEdge === edge)
                    continue;
                const isStart = adjacentEdge.startVertex === this;
                const param = isStart ? adjacentEdge.curve.domain.min : adjacentEdge.curve.domain.max;
                // Check G1 continuity (tangent)
                const tangent1 = newCurve.derivative(edge.startVertex === this ? newCurve.domain.min : newCurve.domain.max, 1).normalize();
                const tangent2 = adjacentEdge.curve.derivative(param, 1).normalize();
                const tangentDeviation = Math.abs(Math.abs(tangent1.dot(tangent2)) - 1);
                if (tangentDeviation > this._tolerance) {
                    return false;
                }
            }
            return true;
        } catch (error) {
            return false;
        }
    }
    // Helper method to update dependent geometry after position change
    updateDependentGeometry() {
        // Update any faces that contain this vertex
        const affectedFaces = new Set();
        for (const edge of this._edges) {
            for (const face of edge.faces) {
                affectedFaces.add(face);
            }
        }
        // Notify faces of geometry update
        for (const face of affectedFaces) {
            if (face.updateGeometry) {
                face.updateGeometry();
            }
        }
    }
    calculateOptimalMergePosition(otherVertex) {
        // Get all unique tangent directions at both vertices
        const tangents1 = this.getUniqueTangentDirections();
        const tangents2 = otherVertex.getUniqueTangentDirections();
        // If either vertex has no tangents, use the other vertex's position
        if (tangents1.length === 0)
            return otherVertex.position.clone();
        if (tangents2.length === 0)
            return this.position.clone();
        // Calculate weights based on edge curvatures
        const weight1 = this.calculateVertexWeight();
        const weight2 = otherVertex.calculateVertexWeight();
        // Calculate weighted position
        const totalWeight = weight1 + weight2;
        if (totalWeight < Number.EPSILON) {
            return this.position.clone();
        }
        return this.position.multiply(weight1).add(otherVertex.position.multiply(weight2)).multiply(1 / totalWeight);
    }
    getUniqueTangentDirections() {
        const tangents = new Set();
        for (const edge of this.edges) {
            const tangent = this.getEdgeTangent(edge);
            // Store normalized tangent vector as string for comparison
            tangents.add(`${ tangent.x },${ tangent.y },${ tangent.z }`);
        }
        return Array.from(tangents);
    }
    calculateVertexWeight() {
        let weight = 0;
        for (const edge of this.edges) {
            try {
                // Use edge curvature at vertex to influence weight
                const param = edge.startVertex === this ? edge.curve.domain.min : edge.curve.domain.max;
                const derivative2 = edge.curve.derivative(param, 2);
                weight += derivative2.length();
            } catch (error) {
                // If curvature calculation fails, use default weight
                weight += 1;
            }
        }
        return weight;
    }
    findMergeableEdgePairs(otherVertex) {
        const pairs = [];
        for (const edge1 of this.edges) {
            for (const edge2 of otherVertex.edges) {
                if (this.areEdgesMergeable(edge1, edge2)) {
                    pairs.push([
                        edge1,
                        edge2
                    ]);
                }
            }
        }
        return pairs;
    }
    areEdgesMergeable(edge1, edge2) {
        // Check if edges share same curve type and are geometrically continuous
        if (edge1.curve.constructor !== edge2.curve.constructor) {
            return false;
        }
        // Check tangent continuity
        const tangent1 = this.getEdgeTangent(edge1);
        const tangent2 = this.getEdgeTangent(edge2);
        const dotProduct = Math.abs(tangent1.dot(tangent2));
        return Math.abs(dotProduct - 1) < this.tolerance;
    }
    mergeEdges(edge1, edge2) {
        // Input validation
        if (!(edge1 instanceof TopologyEdge) || !(edge2 instanceof TopologyEdge)) {
            throw new Error('Both parameters must be TopologyEdge instances');
        }
        // Check if both edges are connected to this vertex
        if (!this.edges.includes(edge1) || !this.edges.includes(edge2)) {
            throw new Error('Both edges must be connected to this vertex');
        }
        // Check if edges can be merged (share this vertex and have compatible curves)
        const isStart1 = edge1.startVertex === this;
        const isStart2 = edge2.startVertex === this;
        try {
            // Create new merged edge
            const mergedEdge = new TopologyEdge();
            // Set up curve orientation and vertex connections based on edge configurations
            if (isStart1 && isStart2) {
                // Both edges start at this vertex - reverse second edge
                edge2.reverse();
                mergedEdge.startVertex = edge2.endVertex;
                mergedEdge.endVertex = edge1.endVertex;
                mergedEdge.curve = edge2.curve.join(edge1.curve);
            } else if (!isStart1 && !isStart2) {
                // Both edges end at this vertex - reverse second edge
                edge2.reverse();
                mergedEdge.startVertex = edge1.startVertex;
                mergedEdge.endVertex = edge2.startVertex;
                mergedEdge.curve = edge1.curve.join(edge2.curve);
            } else if (isStart1 && !isStart2) {
                // First edge starts, second ends at this vertex
                mergedEdge.startVertex = edge2.startVertex;
                mergedEdge.endVertex = edge1.endVertex;
                mergedEdge.curve = edge2.curve.join(edge1.curve);
            } else {
                // First edge ends, second starts at this vertex
                mergedEdge.startVertex = edge1.startVertex;
                mergedEdge.endVertex = edge2.endVertex;
                mergedEdge.curve = edge1.curve.join(edge2.curve);
            }
            // Update vertex-edge references
            mergedEdge.startVertex.edges = mergedEdge.startVertex.edges.filter(e => e !== edge1 && e !== edge2);
            mergedEdge.startVertex.edges.push(mergedEdge);
            mergedEdge.endVertex.edges = mergedEdge.endVertex.edges.filter(e => e !== edge1 && e !== edge2);
            mergedEdge.endVertex.edges.push(mergedEdge);
            // Merge face references
            mergedEdge.faces = Array.from(new Set([
                ...edge1.faces,
                ...edge2.faces
            ]));
            // Update face references to the merged edge
            mergedEdge.faces.forEach(face => {
                const index1 = face.edges.indexOf(edge1);
                const index2 = face.edges.indexOf(edge2);
                if (index1 !== -1)
                    face.edges[index1] = mergedEdge;
                if (index2 !== -1)
                    face.edges.splice(index2, 1);
            });
            // Remove this vertex from topology
            this.edges = [];
            this.position = null;
            // Validate merged edge geometry and topology
            this.validateMergedEdge(mergedEdge);
            // Invalidate original edges
            edge1.curve = null;
            edge1.startVertex = null;
            edge1.endVertex = null;
            edge1.faces = [];
            edge2.curve = null;
            edge2.startVertex = null;
            edge2.endVertex = null;
            edge2.faces = [];
            return mergedEdge;
        } catch (error) {
            throw new Error(`Edge merge operation failed: ${ error.message }`);
        }
    }
    // ... other methods ...
    adjustEdgeGeometry(edge) {
        try {
            if (!(edge instanceof TopologyEdge)) {
                throw new Error('Parameter must be a TopologyEdge instance');
            }
            // Determine if this vertex is the start or end of the edge
            const isStart = edge.startVertex === this;
            const parameter = isStart ? edge.curve.domain.min : edge.curve.domain.max;
            // Handle NURBS curves specifically
            if (edge.curve instanceof NurbsCurve) {
                // Create new control points array
                const controlPoints = [...edge.curve.controlPoints];
                // Update the control point at the appropriate end
                const index = isStart ? 0 : controlPoints.length - 1;
                const weight = controlPoints[index].weight();
                // Create new control point at updated position
                controlPoints[index] = new ControlPoint(this.position.x, this.position.y, this.position.z, weight);
                // Create new curve with adjusted control point
                const newCurve = new NurbsCurve(controlPoints, edge.curve.knotVector, edge.curve.degree);
                // Verify the new curve maintains continuity
                if (this.validateCurveContinuity(newCurve, edge)) {
                    edge.curve = newCurve;
                    return true;
                }
            }
            // If adjustment fails or curve type is unsupported
            return false;
        } catch (error) {
            throw new Error(`Failed to adjust edge geometry: ${ error.message }`);
        }
    }
    sortEdgesByParameter() {
        this.edges.sort((a, b) => {
            const paramA = a.startVertex === this ? a.curve.domain.min : a.curve.domain.max;
            const paramB = b.startVertex === this ? b.curve.domain.min : b.curve.domain.max;
            return paramA - paramB;
        });
    }
    validateMergedEdge(edge) {
        // Verify edge has valid curve
        if (!edge.curve || !(edge.curve instanceof NurbsCurve)) {
            throw new Error('Invalid curve in merged edge');
        }
        // Verify vertex connections
        if (!edge.startVertex || !edge.endVertex) {
            throw new Error('Invalid vertex connections in merged edge');
        }
        // Verify edge is in vertex edge lists
        if (!edge.startVertex.edges.includes(edge) || !edge.endVertex.edges.includes(edge)) {
            throw new Error('Edge not properly referenced in vertex lists');
        }
        // Verify face references
        if (!edge.faces.every(face => face.edges.includes(edge))) {
            throw new Error('Invalid face references in merged edge');
        }
        // Verify curve geometry matches vertex positions
        const startPoint = edge.curve.evaluate(edge.curve.domain.min);
        const endPoint = edge.curve.evaluate(edge.curve.domain.max);
        if (startPoint.subtract(edge.startVertex.position).length() > this.tolerance || endPoint.subtract(edge.endVertex.position).length() > this.tolerance) {
            throw new Error('Edge curve geometry does not match vertex positions');
        }
        // Verify tangent continuity at join point
        const midParam = (edge.curve.domain.min + edge.curve.domain.max) / 2;
        const tangent = edge.curve.derivative(midParam, 1);
        if (tangent.length() < Number.EPSILON) {
            throw new Error('Invalid tangent at curve join point');
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
                throw new Error(`Trim curves must form a closed loop: Gap detected between curves ${ i } and ${ (i + 1) % trimCurves.length }`);
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
                throw new Error(`Trim curves must form a closed loop: Gap detected between curves ${ i } and ${ (i + 1) % trimCurves.length }`);
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
                        throw new Error(`Failed to check curve intersections: ${ error.message }`);
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
        try {
            // Input validation
            if (!trimLoop || !Array.isArray(trimLoop.curves)) {
                throw new Error('Invalid trim loop structure');
            }
            if (trimLoop.curves.length === 0) {
                throw new Error('Empty trim loop');
            }
            // Calculate loop orientation
            let signedArea = 0;
            const curves = trimLoop.curves;
            // Calculate signed area in parameter space
            for (let i = 0; i < curves.length; i++) {
                const curve = curves[i];
                const nextCurve = curves[(i + 1) % curves.length];
                // Sample points at start and end of each curve
                const p1 = curve.evaluate(curve.domain.min);
                const p2 = curve.evaluate(curve.domain.max);
                const p3 = nextCurve.evaluate(nextCurve.domain.min);
                // Calculate contribution to signed area using cross product
                signedArea += (p2.x - p1.x) * (p2.y + p1.y) / 2;
                signedArea += (p3.x - p2.x) * (p3.y + p2.y) / 2;
            }
            // Check if orientation matches the expected orientation
            const currentOrientation = signedArea > 0;
            const expectedOrientation = !trimLoop.isOuter;
            // CCW for holes, CW for outer boundary
            // If orientation is incorrect, reverse all curves in the loop
            if (currentOrientation !== expectedOrientation) {
                for (let i = 0; i < trimLoop.curves.length; i++) {
                    trimLoop.curves[i] = trimLoop.curves[i].reverse();
                }
                // Reverse the order of curves in the loop
                trimLoop.curves.reverse();
            }
            // Verify curve connectivity after potential reversal
            for (let i = 0; i < curves.length; i++) {
                const currentCurve = curves[i];
                const nextCurve = curves[(i + 1) % curves.length];
                // Get end points
                const currentEnd = currentCurve.evaluate(currentCurve.domain.max);
                const nextStart = nextCurve.evaluate(nextCurve.domain.min);
                // Check connectivity within tolerance
                if (currentEnd.subtract(nextStart).length() > this.tolerance) {
                    throw new Error(`Gap detected between trim curves at index ${ i }`);
                }
            }
            return {
                isValid: true,
                orientation: currentOrientation ? 'CCW' : 'CW',
                wasReversed: currentOrientation !== expectedOrientation
            };
        } catch (error) {
            throw new Error(`Trim curve orientation validation failed: ${ error.message }`);
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
    subdivideCell(i, j, gridPoints, gridParams) {
        try {
            // Validate input parameters
            if (!Number.isInteger(i) || !Number.isInteger(j)) {
                throw new Error('Grid indices must be integers');
            }
            if (!Array.isArray(gridPoints) || !Array.isArray(gridParams)) {
                throw new Error('Grid arrays must be valid arrays');
            }
            // Get corner points and parameters of the cell
            const p00 = gridPoints[i][j];
            const p10 = gridPoints[i + 1][j];
            const p01 = gridPoints[i][j + 1];
            const p11 = gridPoints[i + 1][j + 1];
            const params00 = gridParams[i][j];
            const params10 = gridParams[i + 1][j];
            const params01 = gridParams[i][j + 1];
            const params11 = gridParams[i + 1][j + 1];
            // Calculate midpoints parameters
            const midU = {
                u: (params00.u + params10.u) / 2,
                v: (params00.v + params10.v) / 2
            };
            const midV = {
                u: (params00.u + params01.u) / 2,
                v: (params00.v + params01.v) / 2
            };
            const midUV = {
                u: (params00.u + params11.u) / 2,
                v: (params00.v + params11.v) / 2
            };
            // Evaluate surface at midpoints
            const pMidU = this.isPointInside(midU.u, midU.v) ? this.baseSurface.evaluate(midU.u, midU.v) : null;
            const pMidV = this.isPointInside(midV.u, midV.v) ? this.baseSurface.evaluate(midV.u, midV.v) : null;
            const pMidUV = this.isPointInside(midUV.u, midUV.v) ? this.baseSurface.evaluate(midUV.u, midUV.v) : null;
            // Insert new rows and columns in the grids
            gridPoints.splice(i + 1, 0, new Array(gridPoints[i].length).fill(null));
            gridParams.splice(i + 1, 0, new Array(gridParams[i].length).fill(null));
            // Insert new columns
            for (let row = 0; row < gridPoints.length; row++) {
                gridPoints[row].splice(j + 1, 0, null);
                gridParams[row].splice(j + 1, 0, null);
            }
            // Update grid points with new midpoints
            if (pMidU) {
                gridPoints[i + 1][j] = pMidU;
                gridParams[i + 1][j] = midU;
            }
            if (pMidV) {
                gridPoints[i][j + 1] = pMidV;
                gridParams[i][j + 1] = midV;
            }
            if (pMidUV) {
                gridPoints[i + 1][j + 1] = pMidUV;
                gridParams[i + 1][j + 1] = midUV;
            }
            if (p10) {
                gridPoints[i + 2][j] = p10;
                gridParams[i + 2][j] = params10;
            }
            if (p01) {
                gridPoints[i][j + 2] = p01;
                gridParams[i][j + 2] = params01;
            }
            if (p11) {
                gridPoints[i + 2][j + 2] = p11;
                gridParams[i + 2][j + 2] = params11;
            }
            // Fill in any remaining boundary points if needed
            if (p10 && p11) {
                const midUR = {
                    u: (params10.u + params11.u) / 2,
                    v: (params10.v + params11.v) / 2
                };
                if (this.isPointInside(midUR.u, midUR.v)) {
                    gridPoints[i + 2][j + 1] = this.baseSurface.evaluate(midUR.u, midUR.v);
                    gridParams[i + 2][j + 1] = midUR;
                }
            }
            if (p01 && p11) {
                const midLR = {
                    u: (params01.u + params11.u) / 2,
                    v: (params01.v + params11.v) / 2
                };
                if (this.isPointInside(midLR.u, midLR.v)) {
                    gridPoints[i + 1][j + 2] = this.baseSurface.evaluate(midLR.u, midLR.v);
                    gridParams[i + 1][j + 2] = midLR;
                }
            }
        } catch (error) {
            throw new Error(`Cell subdivision failed: ${ error.message }`);
        }
    }
    // ...existing methods...
    clone() {
        try {
            const clonedSurface = new TrimmedSurface();
            clonedSurface.baseSurface = this.baseSurface.clone();
            clonedSurface.trimCurves = this.trimCurves.map(loop => ({
                curves: loop.curves.map(curve => curve.clone()),
                isOuter: loop.isOuter
            }));
            return clonedSurface;
        } catch (error) {
            throw new Error(`Failed to clone trimmed surface: ${ error.message }`);
        }
    }
    compute() {
        // Input validation
        if (!this.baseSurface || !(this.baseSurface instanceof NurbsSurface)) {
            throw new Error('Base surface must be set before computing the trimmed surface');
        }
        // Initialize storage for vertices, edges, and faces
        const vertices = [];
        const edges = [];
        const faces = [];
        // Iterate through trim curves and compute points
        this.trimCurves.forEach(trimLoop => {
            const loopVertices = [];
            const loopEdges = [];
            trimLoop.curves.forEach(curve => {
                const tessellation = curve.tessellate();
                // Assume tessellate is a method that gives us points
                tessellation.points.forEach(point => {
                    const vertex = new TopologyVertex({ position: point });
                    loopVertices.push(vertex);
                    vertices.push(vertex);
                });
                // Create edges from the tessellated points
                for (let i = 0; i < tessellation.points.length - 1; i++) {
                    const edge = new TopologyEdge();
                    edge.startVertex = loopVertices[i];
                    edge.endVertex = loopVertices[i + 1];
                    edge.curve = curve;
                    loopEdges.push(edge);
                    edges.push(edge);
                }
            });
            // Create a face from the edges
            const face = new TopologyFace();
            loopEdges.forEach(edge => {
                face.addBound(edge);
            });
            faces.push(face);
        });
        // Return constructed topology
        return {
            vertices: Object.freeze(vertices),
            edges: Object.freeze(edges),
            faces: Object.freeze(faces)
        };
    }
}
// ... other methods
class Vector3D {
    constructor(x = 0, y = 0, z = 0) {
        // Input validation
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
            throw new Error('Coordinates must be finite numbers');
        }
        // Use Object.defineProperties for immutable coordinates
        Object.defineProperties(this, {
            'x': {
                value: x === 0 ? 0 : x,
                // Protect against -0
                writable: false,
                enumerable: true,
                configurable: false
            },
            'y': {
                value: y === 0 ? 0 : y,
                writable: false,
                enumerable: true,
                configurable: false
            },
            'z': {
                value: z === 0 ? 0 : z,
                writable: false,
                enumerable: true,
                configurable: false
            }
        });
        // Freeze the instance to make it fully immutable
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
    // ... other methods ...
    project(onto) {
        // Input validation
        if (!(onto instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        // Calculate dot product
        const dotProduct = this.dot(onto);
        const ontoLengthSquared = onto.lengthSquared();
        if (Math.abs(ontoLengthSquared) < Number.EPSILON) {
            throw new Error('Cannot project onto zero vector');
        }
        const scale = dotProduct / ontoLengthSquared;
        // Return projected vector
        return new Vector3D(onto.x * scale, onto.y * scale, onto.z * scale);
    }
    reject(from) {
        // Input validation
        if (!(from instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        // Check if 'from' vector is zero vector
        if (Math.abs(from.x) < Number.EPSILON && Math.abs(from.y) < Number.EPSILON && Math.abs(from.z) < Number.EPSILON) {
            throw new Error('Cannot reject from zero vector');
        }
        try {
            // Calculate projection first
            const projection = this.project(from);
            // Rejection is the difference between original vector and projection
            return this.subtract(projection);
        } catch (error) {
            throw new Error(`Vector rejection failed: ${ error.message }`);
        }
    }
    rotate(axis, angleRadians) {
        if (!(axis instanceof Vector3D)) {
            throw new Error('Axis must be a Vector3D instance');
        }
        if (typeof angleRadians !== 'number' || !Number.isFinite(angleRadians)) {
            throw new Error('Angle must be a finite number');
        }
        const cosTheta = Math.cos(angleRadians);
        const sinTheta = Math.sin(angleRadians);
        const {
            x: ax,
            y: ay,
            z: az
        } = axis.normalize();
        const rotationMatrix = [
            [
                cosTheta + ax * ax * (1 - cosTheta),
                ax * ay * (1 - cosTheta) - az * sinTheta,
                ax * az * (1 - cosTheta) + ay * sinTheta
            ],
            [
                ay * ax * (1 - cosTheta) + az * sinTheta,
                cosTheta + ay * ay * (1 - cosTheta),
                ay * az * (1 - cosTheta) - ax * sinTheta
            ],
            [
                az * ax * (1 - cosTheta) - ay * sinTheta,
                az * ay * (1 - cosTheta) + ax * sinTheta,
                cosTheta + az * az * (1 - cosTheta)
            ]
        ];
        const newX = rotationMatrix[0][0] * this.x + rotationMatrix[0][1] * this.y + rotationMatrix[0][2] * this.z;
        const newY = rotationMatrix[1][0] * this.x + rotationMatrix[1][1] * this.y + rotationMatrix[1][2] * this.z;
        const newZ = rotationMatrix[2][0] * this.x + rotationMatrix[2][1] * this.y + rotationMatrix[2][2] * this.z;
        return new Vector3D(newX, newY, newZ);
    }
    // ... other methods ...
    length() {
        const lengthSquared = this.x * this.x + this.y * this.y + this.z * this.z;
        if (lengthSquared < 0) {
            throw new Error('Invalid length calculation: negative length squared');
        }
        return Math.sqrt(lengthSquared);
    }
    lengthSquared() {
        try {
            // Calculate sum of squares for each component
            const squaredLength = this.x * this.x + this.y * this.y + this.z * this.z;
            // Check for numerical overflow/underflow
            if (!Number.isFinite(squaredLength)) {
                throw new Error('Length calculation resulted in non-finite value');
            }
            // Return the squared length
            return squaredLength === 0 ? 0 : squaredLength;
        } // Protect against -0
        catch (error) {
            throw new Error(`Failed to calculate squared length: ${ error.message }`);
        }
    }
    distanceTo(vector) {
        // Input validation
        if (!(vector instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        try {
            // Calculate distance using subtraction and length
            const diff = this.subtract(vector);
            return diff.length();
        } catch (error) {
            throw new Error(`Distance calculation failed: ${ error.message }`);
        }
    }
    equals(vector, tolerance = defaultTolerance) {
        // Input validation
        if (!(vector instanceof Vector3D)) {
            throw new Error('Parameter must be a Vector3D instance');
        }
        if (typeof tolerance !== 'number' || tolerance < 0) {
            throw new Error('Tolerance must be a non-negative number');
        }
        try {
            // Compare each component within tolerance
            const dx = Math.abs(this.x - vector.x);
            const dy = Math.abs(this.y - vector.y);
            const dz = Math.abs(this.z - vector.z);
            // Check if differences are within tolerance
            return dx <= tolerance && dy <= tolerance && dz <= tolerance;
        } catch (error) {
            throw new Error(`Vector comparison failed: ${ error.message }`);
        }
    }
    clone() {
        // Simple clone method to create a new Vector3D instance with the same coordinates
        try {
            return new Vector3D(this.x === 0 ? 0 : this.x, // Protect against -0
            this.y === 0 ? 0 : this.y, this.z === 0 ? 0 : this.z);
        } catch (error) {
            throw new Error(`Failed to clone vector: ${ error.message }`);
        }
    }
}
// Helper function to calculate factorial
function factorial(n) {
    if (n === 0 || n === 1)
        return 1;
    return n * factorial(n - 1);
}
function test() {
}