// Basic data structures
class Point {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    isEqual(otherPoint) {
        if (!(otherPoint instanceof Point)) {
            throw new Error('Argument must be an instance of Point.');
        }
        return this.x === otherPoint.x && this.y === otherPoint.y && this.z === otherPoint.z;
    }
    distance(otherPoint) {
        if (!(otherPoint instanceof Point)) {
            throw new Error('Argument must be an instance of Point.');
        }
        const dx = this.x - otherPoint.x;
        const dy = this.y - otherPoint.y;
        const dz = this.z - otherPoint.z;
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
    translate(dx, dy, dz) {
        this.x += dx;
        this.y += dy;
        this.z += dz;
        return this;
    }
    // Return the updated point for method chaining
    scale(scalar) {
        this.x *= scalar;
        this.y *= scalar;
        this.z *= scalar;
        return this;
    }
    // Return the updated point for method chaining
    toArray() {
        return [
            this.x,
            this.y,
            this.z
        ];
    }
    static fromArray(array) {
        if (!Array.isArray(array) || array.length !== 3) {
            throw new Error('Input must be an array with three elements.');
        }
        return new Point(array[0], array[1], array[2]);
    }
}
class Triangle {
    constructor(points) {
        if (points.length !== 3) {
            throw new Error("A triangle must have exactly 3 points.");
        }
        this.points = points;
        this.plane = this.computePlaneEquation();
    }

    computePlaneEquation() {
        const [A, B, C] = this.points;
        const normal = this.crossProduct(this.subtract(B, A), this.subtract(C, A));
        const D = -this.dotProduct(normal, A);
        return { A: normal.x, B: normal.y, C: normal.z, D: D };
    }

    areCoPlanar(otherTriangle) {
        const { A, B, C, D } = this.plane;
        for (let point of otherTriangle.points) {
            let distance = A * point.x + B * point.y + C * point.z + D;
            if (Math.abs(distance) > 1e-6) {
                return false;
            }
        }
        return true;
    }

    intersectPlanes(otherPlane) {
        const direction = this.crossProduct(
            { x: this.plane.A, y: this.plane.B, z: this.plane.C },
            { x: otherPlane.A, y: otherPlane.B, z: otherPlane.C }
        );

        if (this.isZeroVector(direction)) {
            return null; // Planes are parallel
        }

        // Find a point on the line of intersection
        const pointOnLine = this.solvePlanesIntersectionPoint(this.plane, otherPlane);
        return { point: pointOnLine, direction: direction };
    }

    intersectLineWithTriangle(line) {
        const intersectionPoints = [];
        for (let i = 0; i < 3; i++) {
            const p1 = this.points[i];
            const p2 = this.points[(i + 1) % 3];
            const intersection = this.intersectLineWithSegment(line, p1, p2);
            if (intersection) {
                intersectionPoints.push(intersection);
            }
        }
        return intersectionPoints;
    }

    findIntersection(otherTriangle) {
        if (this.areCoPlanar(otherTriangle)) {
            // Handle co-planar intersection (2D polygon intersection)
            return this.intersectCoPlanarTriangles(otherTriangle);
        }

        const line = this.intersectPlanes(otherTriangle.plane);
        if (!line) {
            return []; // Triangles are parallel and do not intersect
        }

        const points1 = this.intersectLineWithTriangle(line);
        const points2 = otherTriangle.intersectLineWithTriangle(line);
        return this.findOverlapSegment(points1, points2, line);
    }

    findOverlapSegment(points1, points2, line) {
        if (points1.length === 0 || points2.length === 0) {
            return [];
        }

        let [t1Min, t1Max] = this.getParametricRange(points1, line);
        let [t2Min, t2Max] = this.getParametricRange(points2, line);

        let overlapMin = Math.max(t1Min, t2Min);
        let overlapMax = Math.min(t1Max, t2Max);

        if (overlapMin <= overlapMax) {
            return [this.evaluateLine(line, overlapMin), this.evaluateLine(line, overlapMax)];
        }

        return [];
    }

    // Helper Methods
    subtract(p1, p2) {
        return { x: p1.x - p2.x, y: p1.y - p2.y, z: p1.z - p2.z };
    }

    dotProduct(v1, v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    crossProduct(v1, v2) {
        return {
            x: v1.y * v2.z - v1.z * v2.y,
            y: v1.z * v2.x - v1.x * v2.z,
            z: v1.x * v2.y - v1.y * v2.x
        };
    }

    isZeroVector(v) {
        return Math.abs(v.x) < 1e-6 && Math.abs(v.y) < 1e-6 && Math.abs(v.z) < 1e-6;
    }

    solvePlanesIntersectionPoint(plane1, plane2) {
        // Solve the linear system to find a point on the line of intersection
        // This implementation is simplified and may need a robust solver for general use
        return { x: 0, y: 0, z: 0 }; // Placeholder
    }

    intersectLineWithSegment(line, p1, p2) {
        // Calculate intersection of a line with a segment
        return null; // Placeholder
    }

    getParametricRange(points, line) {
        // Calculate the parametric range for the line segment
        return [0, 1]; // Placeholder
    }

    evaluateLine(line, t) {
        return {
            x: line.point.x + t * line.direction.x,
            y: line.point.y + t * line.direction.y,
            z: line.point.z + t * line.direction.z
        };
    }
}
class Edge {
    constructor(points) {
        if (points.length < 2) {
            throw new Error('An edge must have at least two points.');
        }
        this.points = points;
    }
    isEqual(otherEdge) {
        if (!(otherEdge instanceof Edge)) {
            throw new Error('Argument must be an instance of Edge.');
        }
        const thisPointSet = new Set(this.points.map(point => `${ point.x },${ point.y },${ point.z }`));
        const otherPointSet = new Set(otherEdge.points.map(point => `${ point.x },${ point.y },${ point.z }`));
        const hasSameLength = thisPointSet.size === otherPointSet.size;
        const hasMatchingPoints = [...thisPointSet].every(point => otherPointSet.has(point));
        // Check for both direct and reverse equality
        const thisPointsReversed = this.points.slice().reverse();
        const otherPointsReversedSet = new Set(this.points.map(point => `${ otherEdge.points[otherEdge.points.length - 1 - this.points.indexOf(point)].x },${ otherEdge.points[otherEdge.points.length - 1 - this.points.indexOf(point)].y },${ otherEdge.points[otherEdge.points.length - 1 - this.points.indexOf(point)].z }`));
        const hasMatchingPointsReversed = [...thisPointSet].every(point => otherPointsReversedSet.has(point));
        return hasSameLength && (hasMatchingPoints || hasMatchingPointsReversed);
    }
}
// Supports polylines
class Face {
    constructor(mesh, edges = []) {
        this.mesh = mesh;
        // Array of Triangle objects
        this.edges = edges;
    }
}
// Array of Edge objects
class Solid {
    constructor(faces = [], edges = []) {
        this.faces = faces;
        // Array of Face objects
        this.edges = edges;
    }
}
// Shared edges for easier Boolean operations
class BREP {
    constructor() {
        this.points = [];
        this.triangles = [];
        this.edges = [];
        this.faces = [];
        this.solids = [];
    }
    addPoint(x, y, z) {
        const point = new Point(x, y, z);
        this.points.push(point);
        return point;
    }
    addTriangle(p1, p2, p3) {
        const triangle = new Triangle(p1, p2, p3);
        this.triangles.push(triangle);
        return triangle;
    }
    addEdge(points) {
        const edge = new Edge(points);
        this.edges.push(edge);
        return edge;
    }
    addFace(mesh, edges = []) {
        const face = new Face(mesh, edges);
        this.faces.push(face);
        return face;
    }
    addSolid(faces, edges = []) {
        const solid = new Solid(faces, edges);
        this.solids.push(solid);
        return solid;
    }
    validateMesh() {
    }
    performBooleanOperation(solidA, solidB, operationType) {
    }
    extrudeFace(face, depth) {
        const newFaces = [];
        const newEdges = [];
        face.edges.forEach(edge => {
            const extrudedPoints = edge.points.map(point => {
                const extrudedPoint = new Point(point.x, point.y, point.z + depth);
                this.points.push(extrudedPoint);
                return extrudedPoint;
            });
            const newEdge = new Edge(extrudedPoints);
            newEdges.push(newEdge);
            this.edges.push(newEdge);
        });
        const mesh = [];
        for (let i = 0; i < edge.points.length - 1; i++) {
            const p1 = edge.points[i];
            const p2 = edge.points[i + 1];
            const p3 = newEdges[i].points[0];
            const p4 = newEdges[i + 1].points[0];
            const triangle1 = new Triangle(p1, p2, p3);
            const triangle2 = new Triangle(p2, p4, p3);
            mesh.push(triangle1, triangle2);
        }
        const newFace = new Face(mesh, newEdges);
        newFaces.push(newFace);
        this.faces.push(newFace);
        return newFace;
    }
    createFaceFromEdges(edges) {
        if (edges.length < 3) {
            throw new Error('A face must be defined by at least three edges.');
        }
        const triangles = [];
        const edgeMap = new Map();
        edges.forEach(edge => {
            const start = edge.points[0];
            const end = edge.points[edge.points.length - 1];
            if (!edgeMap.has(start)) {
                edgeMap.set(start, []);
            }
            edgeMap.get(start).push(edge);
            if (!edgeMap.has(end)) {
                edgeMap.set(end, []);
            }
            edgeMap.get(end).push(edge);
        });
        const visitedEdges = new Set();
        const loops = [];
        const findLoop = startPoint => {
            const edgeLoop = [];
            let currentPoint = startPoint;
            let firstEdge = null;
            while (edgeMap.has(currentPoint)) {
                const outgoingEdges = edgeMap.get(currentPoint);
                let foundEdge = false;
                for (const edge of outgoingEdges) {
                    if (!visitedEdges.has(edge)) {
                        visitedEdges.add(edge);
                        edgeLoop.push(edge);
                        currentPoint = edge.points[edge.points[edge.points.length - 1] === currentPoint ? 0 : edge.points.length - 1];
                        foundEdge = true;
                        if (!firstEdge)
                            firstEdge = edge;
                        break;
                    }
                }
                if (!foundEdge || currentPoint === firstEdge.points[0] && edgeLoop.length > 0) {
                    break;
                }
            }
            if (edgeLoop.length > 0 && currentPoint === firstEdge.points[0]) {
                return edgeLoop;
            }
            return null;
        };
        edgeMap.forEach((_, startPoint) => {
            if (!visitedEdges.has(startPoint)) {
                const loopEdges = findLoop(startPoint);
                if (loopEdges) {
                    loops.push(loopEdges);
                }
            }
        });
        if (loops.length === 0) {
            throw new Error('No valid loops could be formed from the provided edges.');
        }
        loops.forEach(loop => {
            const firstPoint = loop[0].points[0];
            for (let i = 0; i < loop.length; i++) {
                const edge = loop[i];
                const points = edge.points;
                for (let j = 1; j < points.length - 1; j++) {
                    const triangle = new Triangle(firstPoint, points[j], points[j + 1]);
                    triangles.push(triangle);
                }
            }
        });
        const newFace = new Face(triangles, edges);
        this.faces.push(newFace);
        return newFace;
    }
    // Create a face from a single loop (e.g., a circular loop)
    createFaceFromSingleLoop(points) {
        const triangles = [];
        const firstPoint = points[0];
        for (let i = 1; i < points.length - 1; i++) {
            const triangle = new Triangle(firstPoint, points[i], points[i + 1]);
            triangles.push(triangle);
        }
        const newFace = new Face(triangles, [new Edge(points)]);
        this.faces.push(newFace);
        return newFace;
    }
    // Triangulate a given loop without including holes
    triangulateOuterLoop(loop, triangles) {
        const points = loop.flatMap(edge => edge.points);
        const firstPoint = points[0];
        // Simple triangulation of the outer boundary
        for (let i = 1; i < points.length - 1; i++) {
            const triangle = new Triangle(firstPoint, points[i], points[i + 1]);
            triangles.push(triangle);
        }
    }
    // Mark the inner loops as holes, and ensure no triangles cover these areas
    excludeInnerLoop(loop, triangles) {
        // This method marks points from the inner loop so that they are not included
        // in the final face triangulation. We ensure that triangles formed do not cross
        // into the hole region defined by this loop.
        // Here we can imagine creating constraints that prevent triangles from covering
        // these areas. Since we are not using an external library, we use simple logic
        // to manually filter out triangles that fall within the inner loop area.
        const holePoints = loop.flatMap(edge => edge.points);
        // Filter out triangles that are entirely inside any of the inner loops
        for (let i = triangles.length - 1; i >= 0; i--) {
            const triangle = triangles[i];
            if (this.isTriangleInHole(triangle, holePoints)) {
                triangles.splice(i, 1);
            }
        }
    }
    // Remove the triangle that falls inside the hole
    // Helper function to determine if a triangle falls within a hole
    isTriangleInHole(triangle, holePoints) {
        const [p1, p2, p3] = triangle.points;
        // Simple containment check: if all three points of the triangle are inside the hole
        return this.isPointInPolygon(p1, holePoints) && this.isPointInPolygon(p2, holePoints) && this.isPointInPolygon(p3, holePoints);
    }
    // Helper function to check if a point lies inside a polygon (using ray-casting algorithm)
    isPointInPolygon(point, polygonPoints) {
        let crossings = 0;
        for (let i = 0; i < polygonPoints.length; i++) {
            const p1 = polygonPoints[i];
            const p2 = polygonPoints[(i + 1) % polygonPoints.length];
            if (p1.y > point.y !== p2.y > point.y && point.x < (p2.x - p1.x) * (point.y - p1.y) / (p2.y - p1.y) + p1.x) {
                crossings++;
            }
        }
        return crossings % 2 !== 0;
    }
    calculateLoopArea(loop) {
        let area = 0;
        const n = loop.length;
        for (let i = 0; i < n; i++) {
            const p1 = loop[i].points[0];
            const p2 = loop[(i + 1) % n].points[0];
            area += (p2.x - p1.x) * (p2.y + p1.y);
        }
        return area / 2;
    }
}