class BREPKernel {
    addEdge(edge) {
        this.edges.push(edge);
        edge.mesh.userData.type = 'edge';
        this.entities.push(edge.mesh);
        this.scene.add(edge.mesh);
    }

    removeEdge(edgeId) {
        const edge = this.findEdge(edgeId);
        if (edge) {
            this.manageEntityRemoval(edge, this.edges, edgeId);
        }
    }

    addFace(face) {
        this.faces.push(face);
        face.mesh.userData.type = 'face';
        this.entities.push(face.mesh);
        this.scene.add(face.mesh);
    }

    removeFace(faceId) {
        const face = this.findFace(faceId);
        if (face) {
            this.manageEntityRemoval(face, this.faces, faceId);
        }
    }

    getEdges() {
        return this.edges;
    }

    getFaces() {
        return this.faces;
    }

    findEdge(edgeId) {
        return this.edges.find((edge) => edge.id === edgeId) || null;
    }

    findFace(faceId) {
        return this.faces.find((face) => face.id === faceId) || null;
    }

    constructor(scene, camera) {
        this.scene = scene;
        this.camera = camera || this.createOrthographicCamera();
        this.edges = [];
        this.faces = [];
        this.vertices = [];
        this.entities = [];
        this.currentCameraType = 'orthographic';
        this.controls = null;
        this.hoveredEntity = null;
        this.initControls();
        this.setupEventListeners();
        this.initRaycaster();
        window.addEventListener('resize', this.onWindowResize.bind(this));
    }

    updateEdge(edgeId, newEdgeData) {
        const edgeIndex = this.edges.findIndex((edge) => edge.id === edgeId);
        if (edgeIndex !== -1) {
            this.edges[edgeIndex] = {
                ...this.edges[edgeIndex],
                ...newEdgeData,
            };
        } else {
            throw new Error(`Edge with ID ${edgeId} not found.`);
        }
    }

    updateFace(faceId, newFaceData) {
        const faceIndex = this.faces.findIndex((face) => face.id === faceId);
        if (faceIndex !== -1) {
            this.faces[faceIndex] = {
                ...this.faces[faceIndex],
                ...newFaceData,
            };
        } else {
            throw new Error(`Face with ID ${faceId} not found.`);
        }
    }

    clearEdges() {
        this.edges.forEach((edge) => this.removeEdge(edge.id));
        this.edges = [];
    }

    clearFaces() {
        this.faces.forEach((face) => this.removeFace(face.id));
        this.faces = [];
    }

    countEdges() {
        return this.edges.length;
    }

    countFaces() {
        return this.faces.length;
    }

    addVertex(vertex) {
        this.vertices.push(vertex);
        this.entities.push(vertex.point);
        this.scene.add(vertex.point);
    }

    getVertices() {
        return this.vertices;
    }

    findVertex(vertexId) {
        return this.vertices.find((vertex) => vertex.id === vertexId) || null;
    }

    toggleCamera() {
        if (this.currentCameraType === 'orthographic') {
            this.currentCameraType = 'perspective';
            this.camera = this.createPerspectiveCamera();
        } else {
            this.currentCameraType = 'orthographic';
            this.camera = this.createOrthographicCamera();
        }
        this.camera.updateProjectionMatrix();
        this.initControls();
    }

    clearVertices() {
        this.vertices.forEach((vertex) => this.scene.remove(vertex.point));
        this.vertices = [];
    }

    countVertices() {
        return this.vertices.length;
    }

    onClick(event) {
        this.mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        this.mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
        this.raycaster.setFromCamera(this.mouse, this.camera);
        const intersects = this.raycaster.intersectObjects(this.entities);
        if (intersects.length > 0) {
            const selectedEntity = intersects[0].object;
            this.handleEntityClick(selectedEntity);
        }
    }

    handleEntityClick(entity) {
        const entityType = entity.userData.type;
        console.log(`Clicked on ${entityType}:`, entity);
        this.highlightEntity(entity, true);
    }

    initControls() {
        if (this.controls) {
            this.controls.dispose();
        }
        this.controls = new THREE.TrackballControls(
            this.camera,
            this.scene.domElement,
        );
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.25;
        this.controls.screenSpacePanning = false;
    }

    update() {
        if (this.controls) {
            this.controls.update();
        }
    }

    setupEventListeners() {
        this.onClickBound = this.onClick.bind(this);
        this.onMouseMoveBound = this.onMouseMove.bind(this);
        window.addEventListener('click', this.onClickBound);
        window.addEventListener('mousemove', this.onMouseMoveBound);
    }

    initRaycaster() {
        this.raycaster = new THREE.Raycaster();
        this.mouse = new THREE.Vector2();
    }

    onMouseMove(event) {
        this.mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        this.mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
        this.checkHover();
    }

    checkHover() {
        this.raycaster.setFromCamera(this.mouse, this.camera);
        const intersects = this.raycaster.intersectObjects(this.entities);
        if (intersects.length > 0) {
            const hoveredEntity = intersects[0].object;
            if (this.hoveredEntity !== hoveredEntity) {
                this.highlightEntity(hoveredEntity, true);
                if (this.hoveredEntity) {
                    this.highlightEntity(this.hoveredEntity, false);
                }
                this.hoveredEntity = hoveredEntity;
            }
        } else if (this.hoveredEntity) {
            this.highlightEntity(this.hoveredEntity, false);
            this.hoveredEntity = null;
        }
    }

    highlightEntity(entity, highlight) {
        if (entity.material) {
            entity.material.emissive.set(highlight ? 0x555555 : 0x000000);
            if (highlight) {
                setTimeout(() => this.highlightEntity(entity, false), 1000);
            }
        }
    }

    clearEntities() {
        this.clearEdges();
        this.clearFaces();
        this.clearVertices();
        this.entities = [];
    }

    dispose() {
        this.clearEntities();
        this.controls.dispose();
        this.clearEventListeners();
    }

    manageEntityRemoval(entity, storageArray, entityId) {
        storageArray = storageArray.filter((item) => item.id !== entityId);
        this.scene.remove(entity.mesh);
        this.entities = this.entities.filter((e) => e.id !== entityId);
        entity.mesh.geometry.dispose();
    }

    onWindowResize() {
        const aspect = window.innerWidth / window.innerHeight;
        if (this.currentCameraType === 'perspective') {
            this.camera.aspect = aspect;
        } else {
            this.camera.left = -aspect;
            this.camera.right = aspect;
        }
        this.camera.updateProjectionMatrix();
    }

    createOrthographicCamera() {
        const aspect = window.innerWidth / window.innerHeight;
        const camera = new THREE.OrthographicCamera(
            -aspect,
            aspect,
            1,
            -1,
            0.1,
            1000,
        );
        camera.position.set(5, 5, 5);
        camera.lookAt(0, 0, 0);
        return camera;
    }

    createPerspectiveCamera() {
        const aspect = window.innerWidth / window.innerHeight;
        const camera = new THREE.PerspectiveCamera(75, aspect, 0.1, 1000);
        camera.position.set(5, 5, 5);
        camera.lookAt(0, 0, 0);
        return camera;
    }

    clearEventListeners() {
        window.removeEventListener('click', this.onClickBound);
        window.removeEventListener('mousemove', this.onMouseMoveBound);
    }
}


class NURBSCurve {
    constructor(degree, controlPoints, knots) {
        this.degree = degree;
        this.controlPoints = controlPoints;
        this.knots = knots;
        this.validate();
    }

    validate() {
        if (this.controlPoints.length < this.degree + 1) {
            throw new Error(
                `Not enough control points for the given degree: ${this.degree}`,
            );
        }

    if (this.knots.length !== this.controlPoints.length + this.degree + 1) {
            throw new Error(
                `Invalid knot vector length: expected ${this.controlPoints.length + this.degree + 1}, got ${this.knots.length}`,
            );
        }

    if (
            !this.knots.every(
                (knot, index) => index === 0 || knot >= this.knots[index - 1],
            )
        ) {
            throw new Error('Knot vector must be non-decreasing.');
        }
    }

    calculateBasisFunction(i, k, u) {
        if (k === 0) {
            return u >= this.knots[i] && u < this.knots[i + 1] ? 1 : 0;
        }
        const denom1 = this.knots[i + k] - this.knots[i];
        const denom2 = this.knots[i + k + 1] - this.knots[i + 1];
        const coeff1 = denom1 !== 0 ? (u - this.knots[i]) / denom1 : 0;
        const coeff2 = denom2 !== 0 ? (this.knots[i + k + 1] - u) / denom2 : 0;
        return (
            coeff1 * this.calculateBasisFunction(i, k - 1, u) +
            coeff2 * this.calculateBasisFunction(i + 1, k - 1, u)
        );
    }

    getPoint(u) {
        let point = new THREE.Vector3();
        const n = this.controlPoints.length - 1;
        for (let i = 0; i <= n; i++) {
            const Bi = this.calculateBasisFunction(i, this.degree, u);
            point.add(this.controlPoints[i].clone().multiplyScalar(Bi));
        }
        return point;
    }

    getPoints(numPoints = 100) {
        const points = [];
        const uStep =
            (this.knots[this.knots.length - 1] - this.knots[this.degree]) /
            numPoints;
        for (let i = 0; i <= numPoints; i++) {
            const u = this.knots[this.degree] + i * uStep;
            points.push(this.getPoint(u));
        }
        return points;
    }

    createMesh(material) {
        const geometry = new THREE.BufferGeometry();
        const points = this.getPoints();
        const positions = new Float32Array(points.length * 3);
        points.forEach((point, index) => {
            positions.set([point.x, point.y, point.z], index * 3);
        });
        geometry.setAttribute(
            'position',
            new THREE.BufferAttribute(positions, 3),
        );
        const line = new THREE.Line(
            geometry,
            material ||
                new THREE.LineBasicMaterial({
                    color: 0x0000ff,
                }),
        );
        return line;
    }

    createLineSegments(material) {
        const geometry = new THREE.BufferGeometry();
        const points = this.getPoints();
        const positions = new Float32Array(points.length * 3);
        points.forEach((point, index) => {
            positions.set([point.x, point.y, point.z], index * 3);
        });
        geometry.setAttribute(
            'position',
            new THREE.BufferAttribute(positions, 3),
        );
        const lineSegments = new THREE.LineSegments(
            geometry,
            material ||
                new THREE.LineDashedMaterial({
                    color: 0x0000ff,
                }),
        );
        return lineSegments;
    }
}
