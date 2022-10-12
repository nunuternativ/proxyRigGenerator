import os
import sys
import time
from collections import OrderedDict

import pymel.core as pm
import maya.cmds as mc
import maya.mel as mel
import maya.OpenMaya as om

# logger
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# QT modules
import Qt
from Qt import wrapInstance, QtGui, QtCore, QtWidgets, QtCompat
# monkey patch for Qt compatibility 
if Qt.__binding__ in ('PySide2', 'PyQt5'):
    Qt.QtWidgets.QHeaderView.setResizeMode = Qt.QtWidgets.QHeaderView.setSectionResizeMode

# custom modules
from rftool.utils.ui import maya_win
from rf_utils.ui import load as loadUi
from nuTools import misc
reload(misc)
from nuTools.util import matrixConstraint as matcon

_VERSION = 2.5
_MODULEDIR = os.path.dirname(sys.modules[__name__].__file__).replace('\\', '/')
_UIFILE = '%s/ui.ui' %_MODULEDIR
_WINDOW_NAME = 'proxyRigGenerator_MainWindow'
_APP_NAME = 'Proxy Rig Generator'

TMP_COMBINE_GEO_NAME = 'tmpCombineGeo'
DEFAULT_CUTOBJ_SUFFIX = '_PrCut'
DEFAULT_GEO_SUFFIX = '_PrGeo'
DEFAULT_JNT_SUFFIX = '_Jnt'
DEFAULT_SEARCHFOR = '_L_'
DEFAULT_REPLACEWITH = '_R_'
DEFAULT_SCALE = 1.0

form_class, base_class = loadUi.loadUiType(_UIFILE)

# --- TO DO
#// - optimize generate cut proxy
#// - new UI
#// - mirror cut cube
#// - snap select cut cube to joints
#// - convert to skin weights
#// - Fix bug in additive mode
#// - Auto length for new cut cube
#// - Change suffix names

def cutMesh(mesh, cutCube, cutCubeFaceIndices=None):
    # if cutCubeFaceIndices not passed, do all faces on the cube
    if cutCubeFaceIndices == None:
        cutCubeFaceIndices = range(cutCube.numFaces())

    meshName = mesh.shortName()
    cubeMDag = cutCube.__apimdagpath__()
    mitPolygon = om.MItMeshPolygon(cubeMDag)
    while not mitPolygon.isDone():
        # for fi in cutCubeFaceIndices:
        fi = mitPolygon.index()
        if fi in cutCubeFaceIndices:
            center = mitPolygon.center(om.MSpace.kWorld)
            
            # get tri points
            tPoints = om.MPointArray()
            mitPolygon.getPoints(tPoints, om.MSpace.kWorld)

            normal = om.MVector()
            mitPolygon.getNormal(normal, om.MSpace.kWorld)
            # invert it
            normal *= -1

            # Get binormal
            binormal = tPoints[1] - tPoints[0]
            binormal.normalize()

            # Get tangent
            tangent = normal ^ binormal
            tangent.normalize()

            # Create matrix
            matrix = pm.dt.TransformationMatrix(
                [binormal.x, binormal.y, binormal.z, 0.0],
                [tangent.x, tangent.y, tangent.z, 0.0],
                [normal.x, normal.y, normal.z, 0.0],
                [center.x, center.y, center.z, 1.0]
            )
            ro = matrix.getRotation()
            ro.setDisplayUnit(unit='degrees')

            # loc = pm.spaceLocator()
            # pm.xform(loc, ws=True, t=[center.x, center.y, center.z])
            # pm.xform(loc, ws=True, ro=[ro.x, ro.y, ro.z])
            
            mc.polyCut(meshName, 
                ch=False,
                pc=[center.x, center.y, center.z],
                ro=[ro.x, ro.y, ro.z])

        mitPolygon.next()

    # mc.select(meshName, r=True)
    # mel.eval('polyCleanupArgList 4 { "0","1","0","0","1","0","0","0","1","1e-05","0","1e-05","0","1e-05","0","-1","0","1" };')
    # mel.eval('polySoftEdge -a 180 -ch 0 %s' %meshName)
    return mesh

class ProxyRigGenerator(form_class, base_class):
    def __init__(self, parent=None):
        super(ProxyRigGenerator, self).__init__(parent)
        self.setupUi(self)

        # --- UI tweaks
        self.setWindowTitle(QtWidgets.QApplication.translate(_WINDOW_NAME, 
            "{0} - v.{1}".format(_APP_NAME, _VERSION), None))

        # Set default UI 
        self.size_doubleSpinBox.setValue(DEFAULT_SCALE)
        self.cubeSuffix_lineEdit.setText(DEFAULT_CUTOBJ_SUFFIX)
        self.geoSuffix_lineEdit.setText(DEFAULT_GEO_SUFFIX)
        self.jntSuffix_lineEdit.setText(DEFAULT_JNT_SUFFIX)
        self.searchFor_lineEdit.setText(DEFAULT_SEARCHFOR)
        self.replaceWith_lineEdit.setText(DEFAULT_REPLACEWITH)
        self.statusBar.showMessage('Ready.')

        # signals
        self.newCube_pushButton.clicked.connect(self.createProxyCubeUI)
        self.snapCube_pushButton.clicked.connect(self.snapProxyCubeUI)
        self.mirrorCube_pushButton.clicked.connect(self.mirrorCubeUI)
        self.loadGeo_pushButton.clicked.connect(self.loadGeo)
        self.skinGeo_pushButton.clicked.connect(self.loadSkinGeo)
        self.rootJnt_pushButton.clicked.connect(self.loadRootJnt)
        self.cut_pushButton.clicked.connect(self.createProxyRigUI)
        self.mirrorCut_pushButton.clicked.connect(self.mirrorCutUI)
        self.generateSkinWeight_pushButton.clicked.connect(self.generateSkinWeightUI)

        # internal vars
        self.geos = []
        self.skinGeo = None

    def getSelPly(self):
        geos = [s for s in misc.getSel(num='inf', selType='transform') if misc.checkIfPly(s)]
        return geos

    def getSelJnt(self):
        jnts = misc.getSel(num='inf', selType='joint')
        return jnts

    def loadGeo(self):
        geos = self.getSelPly()
        self.geos = []
        if geos:
            self.geos = geos
        else:
            self.proxyGeo_lineEdit.setText('')

        geoNames = [str(s.nodeName()) for s in self.geos]
        names = ', '.join(geoNames)
        self.proxyGeo_lineEdit.setText(names)
        self.proxyGeo_lineEdit.setToolTip('%s geometries.\n%s' %(len(geoNames), '\n'.join(geoNames)))

    def loadSkinGeo(self):
        skinGeos = self.getSelPly()
        if skinGeos:
            if len(skinGeos) > 1:
                logger.warning('Only 1 skin geo allow at a time, choosing first one.')
            self.skinGeo = skinGeos[0]
            self.skinGeo_lineEdit.setText(self.skinGeo.shortName())
        else:
            self.skinGeo = None
            self.skinGeo_lineEdit.setText('')

    def loadRootJnt(self):
        selJnt = misc.getSel(num=1, selType='joint')
        if selJnt:
            self.rootJnt = selJnt
            self.rootJnt_lineEdit.setText(self.rootJnt.shortName())
        else:
            self.rootJnt = None
            self.rootJnt_lineEdit.setText('')

    def createProxyCubeUI(self):
        jnts = self.getSelJnt()
        if not jnts:
            logger.error('Select some joint!')
            return

        # get UI inputs
        scale = float(self.size_doubleSpinBox.value())
        jntSuffix = str(self.jntSuffix_lineEdit.text())
        cubeSuffix = str(self.cubeSuffix_lineEdit.text())
        autoStretch = self.autoLen_checkBox.isChecked()
        jntAxis = str(self.axis_comboBox.currentText())
        if self.base_radioButton.isChecked():
            centered = False
        else:
            centered = True

        pm.undoInfo(openChunk=True)
        self.createProxyCube(jnts, scale, jntSuffix, cubeSuffix, centered, autoStretch, jntAxis)
        pm.undoInfo(closeChunk=True)

    def getConnectMethod(self):
        connectMethod = None
        if self.constraint_radioButton.isChecked():
            connectMethod = 'constraint'
        elif self.parent_radioButton.isChecked():
            connectMethod = 'constraint'
        elif self.matrix_radioButton.isChecked():
            connectMethod = 'matrix'
        return connectMethod

    def createProxyRigUI(self):
        cubes = self.getSelPly()
        if not cubes:
            logger.error('Select cut cubes!')
            return
        if not self.geos:
            logger.error('No geometry found!')
            return
        err_sel = [c for c in cubes if c in self.geos]
        if err_sel:
            logger.error('Cannot cut geo by itself: %s' %err_sel)
            return

        # get UI inputs
        jntSuffix = str(self.jntSuffix_lineEdit.text())
        # cubeSuffix = str(self.cubeSuffix_lineEdit.text())
        geoSuffix = str(self.geoSuffix_lineEdit.text())
        removeNamespace = self.removeNamespace_checkBox.isChecked()
        connectMethod = self.getConnectMethod()

        preSels = pm.selected()
        pm.waitCursor(st=True)
        pm.undoInfo(openChunk=True)

        self.createProxyRig(self.geos, cubes, geoSuffix, jntSuffix, connectMethod, removeNamespace)
        
        pm.waitCursor(st=False)
        pm.select(preSels, r=True)
        pm.undoInfo(closeChunk=True)

    def snapProxyCubeUI(self):
        cutCubes = self.getSelPly()
        jnts = self.getSelJnt()
        if cutCubes and jnts:
            pm.undoInfo(openChunk=True)
            self.snapProxyCubePairs(cutCubes, jnts)
            pm.undoInfo(closeChunk=True)
        elif cutCubes:
            pm.undoInfo(openChunk=True)
            self.snapProxyCubeByName(cutCubes)
            pm.undoInfo(closeChunk=True)
        else:
            logger.error('Select cut cubes to snap to matching joints, or select pairs of cube and joint respectively.')

    def mirrorCubeUI(self):
        cutCubes = self.getSelPly()
        if not cutCubes:
            logger.error('Select cut cubes!')
            return

        searchFor = str(self.searchFor_lineEdit.text())
        replaceWith = str(self.replaceWith_lineEdit.text())

        pm.undoInfo(openChunk=True)
        miCubes = self.mirrorGeos(cutCubes, searchFor, replaceWith)
        pm.undoInfo(closeChunk=True)

        return miCubes

    def mirrorCutUI(self):
        geos = self.getSelPly()
        if not geos:
            logger.error('Select cut geos!')
            return

        suffix = str(self.geoSuffix_lineEdit.text())
        searchFor = str(self.searchFor_lineEdit.text())
        replaceWith = str(self.replaceWith_lineEdit.text())
        connectMethod = self.getConnectMethod()

        pm.undoInfo(openChunk=True)
        self.mirrorCut(geos, suffix, searchFor, replaceWith, connectMethod)
        pm.undoInfo(closeChunk=True)
    
    def mirrorCut(self, geos, suffix, searchFor, replaceWith, connectMethod):
        miGeos = self.mirrorGeos(geos, searchFor, replaceWith)
        for geo in miGeos:
            # delete constraint or any transform, if any
            children = geo.getChildren(type='transform')
            if children:
                pm.delete(children)

            # find matching joint
            parent = self.getMatchingJnt(cube=geo, suffix=suffix)
            if parent:
                geoName = geo.nodeName()
                if connectMethod == 'constraint':
                    pm.parentConstraint(parent, geo, mo=True, n='%s_parCons' %geoName)
                    pm.scaleConstraint(parent, geo, mo=True, n='%s_scaCons' %geoName)
                elif connectMethod == 'parent':
                    pm.parent(geo, parent)
                elif connectMethod == 'matrix':
                    mcons = matcon.MatrixConstraint(parents=[parent], child=geo, mo=True, t=True, r=True, s=True)
                    mcons.doIt()

    def generateSkinWeightUI(self):
        sels = self.getSelPly()
        if not sels:
            logger.error('Select skin cubes!')
            return
        if not self.skinGeo:
            logger.error('No skin geo loaded!')
            return
        if not self.rootJnt:
            logger.error('No root joint loaded!')
            return
        
        cubes, jnts = [], []
        for sel in sels:
            jnt = self.getMatchingJnt(sel)
            if jnt:
                cubes.append(sel)
                jnts.append(jnt)
        if not cubes or not jnts:
            logger.error('Failed getting coresponding joints.')
            return

        additive = self.floodRoot_checkBox.isChecked()

        preSels = pm.selected()
        pm.waitCursor(st=True)
        # pm.undoInfo(openChunk=True)
        self.generateSkinWeight(self.skinGeo, self.rootJnt, cubes, jnts, additive)
        # pm.undoInfo(closeChunk=True)
        pm.select(preSels, r=True)
        pm.waitCursor(st=False)

    def generateSkinWeight(self, geo, rootJnt, cubes, jnts, additive=True):
        self.statusBar.showMessage('Preparing...')
        st = time.time()

        triNodes = []
        geoShp = geo.getShape(ni=True)
        geoMDagPath = geoShp.__apimdagpath__()
        hitPoint = om.MFloatPoint()

        # ------------------------------------
        # get current skinCluster node. If none, create new one
        allSkinJnts = [rootJnt] + jnts
        skc = misc.findRelatedSkinCluster(obj=geo)
        if skc:
            currInfs = pm.skinCluster(skc, q=True, inf =True)
            toAdd = [j for j in allSkinJnts if j not in currInfs]
            pm.skinCluster(geo, e=True, ai=toAdd, lw=True, wt=0.0)
        else:
            skc = pm.skinCluster(geo, allSkinJnts, tsb=True)

        skcName = skc.nodeName()

        # unlock all joints
        allInfs = pm.skinCluster(skc, q=True, inf=True)
        infCount = len(allInfs)
        # misc.setJntLockInfluenceWeights(jnts=allInfs, value=False)

        # flood entire geo skin weights with root jnt
        rootIndex = skc.indexForInfluenceObject(rootJnt)
        numVerts = geo.numVertices()
        weights = OrderedDict()
        vtxWeights = OrderedDict()
        for j in allInfs:
            index = skc.indexForInfluenceObject(j)
            value = 0.0
            if index == rootIndex:
                value = 1.0
            vtxWeights[index] = value

        if not additive:
            for vi in xrange(numVerts):
                weights[vi] = vtxWeights

        # ------------------------------------
        # iterate over each skin cube, 
        # assign weights on vertices inside the cube to coresponding joint
        self.statusBar.showMessage('Generating skin weights - 0%')
        ni = 0
        numCubes = float(len(cubes))
        for cube, jnt in zip(cubes, jnts):
            cubeShp = cube.getShape(ni=True)
            jntIndex = skc.indexForInfluenceObject(jnt)

            # get cube data
            cubeMfnMesh, cubeNormals, cubeMBB, tri_node = self.getCutCubeData(cubeShp)
            cubeAccelParams = cubeMfnMesh.autoUniformGridParams()
            if tri_node:
                triNodes.append(tri_node)

            # get point inside cube BB
            mitVtx = om.MItMeshVertex(geoMDagPath)
            # vidInsideCube = []
            while not mitVtx.isDone():
                pt = mitVtx.position(om.MSpace.kWorld)
                if cubeMBB.contains(pt) == True:
                    fPt = om.MFloatPoint(pt[0], pt[1], pt[2])
                    # intersect each normal on cut cube
                    for normal in cubeNormals:
                        intersect = self.rayIntersect(cubeMfnMesh, hitPoint, fPt, normal, cubeAccelParams)
                        # not intersecting even only 1 normal means the point is outside the cut cube
                        if not intersect:
                            break
                    else:
                        vid = mitVtx.index()
                        newWeights = OrderedDict()
                        for ij in vtxWeights:
                            value = 0.0
                            if ij == jntIndex:
                                value = 1.0
                            newWeights[ij] = value
                        weights[vid] = newWeights
                        # vidInsideCube.append(vid)

                mitVtx.next()

            ni += 1
            percentage = (ni/numCubes) * 100.0
            self.statusBar.showMessage('Generating skin weights - %s%%' %(round(percentage, 1)))

        self.statusBar.showMessage('Setting skin weights - 0%')
        numWeights = len(weights)
        ni = 0
        # print weights
        for vid, values in weights.iteritems():
            for jid, value in values.iteritems():
                wlAttr = '%s.weightList[%s].weights[%s]' %(skcName, vid, jid)
                # print vid, jid, values[jid], wlAttr
                mc.setAttr(wlAttr, value)
            if ni % 100 == 0:
                percentage = (vid/float(numVerts)) * 100.0
                self.statusBar.showMessage('Setting skin weights - %s%%' %(round(percentage, 1)))
            ni += 1

        self.statusBar.showMessage('Cleaning up...')
        if triNodes:
            mc.delete(triNodes)

        logger.info('Time spent to generate: %s' %(time.time() - st))
        self.statusBar.showMessage('Ready.')

    def getCutCubeData(self, cutCube):
        cutCubeName = cutCube.shortName()
        mSel = om.MSelectionList()
        mSel.add(cutCubeName)
        cubeMDagPath = om.MDagPath()
        mSel.getDagPath(0, cubeMDagPath)

        tri_node = None
        faceIt = om.MItMeshPolygon(cubeMDagPath)
        nonPlanarFids = []
        while not faceIt.isDone():
            if not faceIt.isPlanar():
                nonPlanarFids.append(faceIt.index())
            faceIt.next()

        # triangulate non planar faces
        if nonPlanarFids:
            nonPlanarFaces = misc.component_range_merge(geoName=cutCubeName, 
                                                    inputList=nonPlanarFids, 
                                                    componentType='face')
            tri_node = mc.polyTriangulate(nonPlanarFaces, ch=True)[0]

        # iterate thru each face of the cube, get face normal and invert it
        faceIt = om.MItMeshPolygon(cubeMDagPath)
        normals = []

        while not faceIt.isDone():
            normalMVector = om.MVector()
            faceIt.getNormal(normalMVector, om.MSpace.kWorld)

            # convert MVector to MFloat vector
            normalMFVector = om.MFloatVector(normalMVector)
            normals.append(normalMFVector)

            # go next
            faceIt.next()

        # get mfn Mesh
        cubeMfnMesh = om.MFnMesh(cubeMDagPath)

        # get inclusive matrix of the cube full path
        cubeIncMatrix = cubeMDagPath.inclusiveMatrix()

        # get cube BoundingBox
        cubeDagFn = om.MFnDagNode(cubeMDagPath)
        cubeMBB = cubeDagFn.boundingBox()  # get bounding box in object space
        cubeMBB.transformUsing(cubeIncMatrix)  # transform the bounding box to world space

        return cubeMfnMesh, normals, cubeMBB, tri_node

    def createProxyCube(self, jnts, mult, searchSuffix, replaceSuffix, centered=True, autoStretch=True, jntAxis='+y'):
        cubes = []
        numJnt = len(jnts)

        # get aim faces
        fdict = {'+x':4, '-x':5, '+y':1, '-y':3, '+z':0, '-z':2}
        oppositeFids = [(4, 5), (1, 3), (0, 2)]
        aimFid = fdict[jntAxis]
        opFid = None
        for fids in oppositeFids:
            if aimFid in fids:
                index = fids.index(aimFid)
                opFid = fids[int(not(index))]
                break

        # get scale vector
        scaleIndex = 0
        for i, axis in enumerate('xyz'):
            if jntAxis.endswith(axis):
                scaleIndex = i
        scaleValue = [mult, mult, mult]
        scaleValue[scaleIndex] = 1

        doMoveFace = True
        if autoStretch and numJnt == 1:
            children = jnts[0].getChildren(type='joint')
            if len(children) > 1 or not children:
                doMoveFace = False
                jnts.append(jnts[0])
                numJnt = 2
                scaleValue[scaleIndex] = mult
            else:
                doMoveFace = True
                jnts.append(children[0])
                numJnt = 2

        lastJnt = None
        for ji, j in enumerate(jnts):
            if autoStretch:
                # skipping last joint
                if ji < (numJnt - 1):
                    cube = pm.polyCube(ch=False, w=1, h=1, d=1)[0]
                    if not centered:
                        # move pivot to bottom
                        pm.xform(cube, ws=True, piv=[0, -0.5, 0])

                    # snap to joint
                    misc.snapTransform('parent', j, cube, False, True)

                    # scale faces on all axis except for the aim axis
                    pm.scale(cube.f, scaleValue, r=True, os=True)

                    # move the faces
                    if doMoveFace == True:
                        nextJnt = jnts[ji+1]
                        nextJntPos = pm.dt.Point(pm.xform(nextJnt, q=True, ws=True, t=True))
                        face = cube.f[aimFid]
                        faceCenter = pm.dt.Point(face.__apimfn__().center(om.MSpace.kWorld))

                        if centered:
                            fMult = 0.5
                            if ji == (numJnt - 2):
                                fMult = 1.0
                            currJntPos = pm.dt.Point(pm.xform(j, q=True, ws=True, t=True))
                            tr = (((nextJntPos - currJntPos) * fMult) + currJntPos) - faceCenter
                            # print tr, ((nextJntPos - currJntPos) * 0.5)
                            pm.move(face, tr, r=True)
                            # break
                            opface = cube.f[opFid]
                            opFaceCenter = pm.dt.Point(opface.__apimfn__().center(om.MSpace.kWorld))
                            if lastJnt:
                                opPos = pm.dt.Point(pm.xform(lastJnt, q=True, ws=True, t=True))
                                tr = (((opPos - currJntPos) * 0.5) + currJntPos) - opFaceCenter
                            else:
                                tr *= -1
                            pm.move(opface, tr, r=True)

                            lastJnt = j
                        else:
                            tr = nextJntPos - faceCenter
                            pm.move(face, tr, r=True)
                else:
                    break
            else:
                cube = pm.polyCube(ch=False, w=mult, h=mult, d=mult)[0]
                misc.snapTransform('parent', j, cube, False, True)

            newName = j.nodeName().replace(searchSuffix, replaceSuffix)
            cube.rename(newName)
            cubes.append(cube)

        pm.select(cubes)
        return cubes

    def getMatchingJnt(self, cube, suffix=None):
        if not suffix:
            suffix = str(self.cubeSuffix_lineEdit.text())
        jntSuffix = str(self.jntSuffix_lineEdit.text())

        cutCubeName = cube.nodeName()
        # if cutCubeName.endswith(cubeSuffix):
        if suffix in cutCubeName:
            jntName = cutCubeName.replace(suffix, jntSuffix)
            jntName = jntName.split(':')[-1]
            for s in ('*:%s' %jntName, jntName):
                jntLs = pm.ls(s, type='joint')
                if jntLs:
                    return jntLs[0]
            else:
                logger.warning('%s: could not find matching joint name.' %cutCubeName)

    def snapProxyCubePairs(self, cutCubes, jnts):
        for cube, jnt in zip(cutCubes, jnts):
            misc.unlockChannelbox([cube], False)
            misc.snapTransform('parent', jnt, cube, False, True)

    def snapProxyCubeByName(self, cutCubes):
        for cube in cutCubes:
            jnt = self.getMatchingJnt(cube)
            if jnt:
                misc.unlockChannelbox([cube], False)
                misc.snapTransform('parent', jnt, cube, False, True)

    def mirrorGeos(self, geos, searchFor, replaceWith):
        miGeos = []
        for geo in geos:
            miGeo = pm.duplicate(geo)[0]
            misc.unlockChannelbox([miGeo], False)

            # move object
            currTr = geo.getTranslation(space='world')
            # scale vertices
            pm.move(miGeo, [currTr.x*-1, currTr.y, currTr.z], ws=True)
            pm.scale(miGeo.vtx, [-1, 1, 1], ws=True, ocp=True)

            # reverse normal
            pm.polyNormal(miGeo, normalMode=0, userNormalMode=0, ch=False)

            # rename
            newName = geo.nodeName().replace(searchFor, replaceWith)
            miGeo.rename(newName)
            miGeos.append(miGeo)

        pm.select(miGeos, r=True)
        return miGeos

    def rayIntersect(self, fnMesh, hitPoint, raySource, rayDir, accelParams=None):
        # hitPoints = om.MFloatPoint()
        # hitRayParams = om.MFloat()
        # hitFaces = om.MIntArray()
        hit = fnMesh.anyIntersection(raySource,  # raySource 
                rayDir,     # rayDirection 
                None,       # faceIds 
                None,       # triIds 
                True,       # idsSorted 
                om.MSpace.kWorld,  # space
                10000,     # maxParam 
                False,      # testBothDirections 
                accelParams,  # accelParams
                hitPoint,  # hitPoints
                None,  # hitRayParams 
                None,  # hitFaces
                None,  # hitTriangles
                None,  # hitBary1s
                None,  # hitBary2s
                1e-03)  # tolerance

        return hit

    def createProxyRig(self, srcGeos, cutCubes, geoSuffix, jntSuffix, connectMethod='constraint', removeNamespace=True):
        self.statusBar.showMessage('Preparing...')
        st = time.time()

        # duplicate and delete history
        dups = pm.duplicate(srcGeos)
        pm.delete(dups, ch=True)

        # unlock channelbox
        misc.unlockChannelbox(dups, False)

        # unhide
        pm.showHidden(dups)

        # unparent to world
        for dup in dups:
            if dup.getParent():
                pm.parent(dup, w=True)

        # get each geo bounding boxes
        geoBBs = {}
        for geo in dups:
            geoShp = geo.getShape(ni=True)
            geoMDagPath = geoShp.__apimdagpath__()
            # get inclusive matrix of the cube full path
            geoIncMatrix = geoMDagPath.inclusiveMatrix()

            # get BoundingBox
            geoDagFn = om.MFnDagNode(geoMDagPath)
            geoMBB = geoDagFn.boundingBox()  # get bounding box in object space
            geoMBB.transformUsing(geoIncMatrix)  # transform the bounding box to world space
            geoBBs[geo] = geoMBB

        # loop over all cut objects
        proxyGeos = []
        triNodes = []
        outsideFaces = []
        newNames = {}  # {pynode:new_name}
        hitPoint = om.MFloatPoint()
        numCubes = float(len(cutCubes))

        ni = 0
        for cutCube in cutCubes:
            cutCubeName = cutCube.nodeName()
            percentage = (ni/numCubes) * 100.0
            self.statusBar.showMessage('Cutting %s - %s%%' %(cutCubeName, round(percentage, 1)))

            # --- find matching joint for cut cube by name
            cutCubeShp = cutCube.getShape(ni=True)
            parent = self.getMatchingJnt(cutCube)
            
            # --- get cut cube data
            cubeMfnMesh, cubeNormals, cubeMBB, tri_node = self.getCutCubeData(cutCubeShp)
            if tri_node:
                triNodes.append(tri_node)

            # --- prepare the mesh
            # loop over each geo finds if the shape intersects the bounding box of the cut cube
            geoInsides = []
            for geo, geoBB in geoBBs.iteritems():
                if geoBB.intersects(cubeMBB):
                    geoInsides.append(geo)

            if geoInsides:
                # duplicate the geos that has BB intersecting box BB and combine them
                newGeos = pm.duplicate(geoInsides)
                if len(newGeos) > 1:
                    proxyTr = pm.polyUnite(newGeos, ch=False, mergeUVSets=True, name=TMP_COMBINE_GEO_NAME)[0]
                    try:
                        pm.delete(newGeos)  # delete the left over transform
                    except:
                        pass
                else:
                    proxyTr = newGeos[0]

                # get cmbGeo data
                proxyGeoShp = proxyTr.getShape(ni=True)
                proxyMDagPath = proxyGeoShp.__apimdagpath__()
                faceCount = proxyGeoShp.numFaces()
                
                # get faces inside cube BB
                mitVtx = om.MItMeshVertex(proxyMDagPath)
                insideBBIndices = set()
                while not mitVtx.isDone():
                    pt = mitVtx.position(om.MSpace.kWorld)
                    if cubeMBB.contains(pt) == True:
                        con_faces_indices = om.MIntArray()
                        mitVtx.getConnectedFaces(con_faces_indices)
                        for ci in xrange(con_faces_indices.length()):
                            cf_fid = con_faces_indices[ci]
                            insideBBIndices.add(cf_fid)

                    mitVtx.next()

                # invert the face id to get the id of faces outside cubeBB
                outsideBBIndices = [i for i in xrange(faceCount) if i not in insideBBIndices]
                outsideBBFaces = misc.component_range_merge(geoName=proxyGeoShp.shortName(), 
                                                        inputList=outsideBBIndices, 
                                                        componentType='face')
                # delete face outside cube BB
                mc.delete(outsideBBFaces)

                # cut the mesh
                cutGeo = cutMesh(mesh=proxyGeoShp, 
                        cutCube=cutCubeShp, 
                        cutCubeFaceIndices=None)
                cutTr = cutGeo.getParent()
                proxyGeos.append(cutTr)

                # iterate each face of the cut geo and see if face center is outside cut cube
                outsideIndices = set()
                cutGeoMDagPath = cutGeo.__apimdagpath__()
                faceIt = om.MItMeshPolygon(cutGeoMDagPath)
                # create accel param
                cubeAccelParams = cubeMfnMesh.autoUniformGridParams()

                while not faceIt.isDone():
                    # get face center
                    centerMPt = faceIt.center(om.MSpace.kWorld)
                    centerMFPt = om.MFloatPoint(centerMPt[0], centerMPt[1], centerMPt[2])

                    # intersect each normal on cut cube
                    for normal in cubeNormals:
                        intersect = self.rayIntersect(cubeMfnMesh, hitPoint, centerMFPt, normal, cubeAccelParams)
                        # not intersecting even only 1 normal means the point is outside the cut cube
                        if not intersect:
                            outsideIndices.add(faceIt.index())
                            break

                    faceIt.next()
                
                # collect faces to be deleted once we exit the loop
                if outsideIndices:
                    extraFaces = misc.component_range_merge(geoName=cutGeo.shortName(), 
                                                        inputList=outsideIndices, 
                                                        componentType='face')
                    outsideFaces.extend(extraFaces)

                # parent cut geo
                if parent:
                    parentName = parent.nodeName()
                    geoName = parentName.replace(jntSuffix, geoSuffix)
                    if removeNamespace:
                        geoName = geoName.split(':')[-1]
                    newNames[cutTr] = geoName  # collect names for cut geos but don't rename them yet

                    if connectMethod == 'constraint':
                        pm.parentConstraint(parent, cutTr, mo=True, n='%s_parCons' %geoName)
                        pm.scaleConstraint(parent, cutTr, mo=True, n='%s_scaCons' %geoName)
                    elif connectMethod == 'parent':
                        pm.parent(cutTr, parent)
                    elif connectMethod == 'matrix':
                        mcons = matcon.MatrixConstraint(parents=[parent], child=cutTr, mo=True, t=True, r=True, s=True)
                        mcons.doIt()
            ni += 1
            
        # cleanup extra faces, triangulat nodes and the master duplicate geos
        self.statusBar.showMessage('Cleaning up...') 
        mc.delete(triNodes)
        mc.delete(outsideFaces)
        mc.delete([n.shortName() for n in dups])

        # rename cut geo to proper name
        for cutGeo, newName in newNames.iteritems():
            cutGeo.rename(newName)

        logger.info('Time spent to generate: %s' %(time.time() - st))
        self.statusBar.showMessage('Ready.')

def show():
    ''' Launching function '''
    global ProxyRigGeneratorApp
    try:
        ProxyRigGeneratorApp.close()
    except Exception, e:
        pass
    
    ProxyRigGeneratorApp = ProxyRigGenerator(parent=maya_win.getMayaWindow())
    ProxyRigGeneratorApp.show()
    return ProxyRigGeneratorApp

'''
from nuTools.util.proxyRigGenerator import app as prgapp
reload(prgapp)

prgapp.show()
'''