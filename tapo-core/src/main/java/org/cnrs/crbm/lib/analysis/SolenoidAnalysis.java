package org.cnrs.crbm.lib.analysis;

import org.apache.commons.math3.geometry.euclidean.threed.Plane;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.geometry.euclidean.twod.Vector2D;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.module.ProVector;
import org.cnrs.crbm.lib.repeats.module.VectorModule;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by pdoviet on 7/10/2015.
 * This class aims to determine right-handed or left-handed solenoid structure.
 */
public class SolenoidAnalysis {

    public static void main(String[] args) throws Exception {


        // 1dce_A
        SolenoidAnalysis solenoidAnalysis = new SolenoidAnalysis();
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<RowRepeatDB> rowsRB = csvReader.getRepeatDB("data/RDB-dataset.tab");
//
//        List<Row> rows = new ArrayList<Row>();
//        //convert to rows
        for (RowRepeatDB r : rowsRB) {

            if (r.getStrclass().equals("III.3")) {
                try {

                    int startRegion = 0;
                    int endRegion = 0;
                    String region = r.getRegion();
                    if (!region.equals("")) {
                        //System.out.println(r.getEntry() + ":" + region);
                        if (region.startsWith("-")) {
                            region = region.substring(1, region.length());
                            startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                            endRegion = Integer.parseInt(region.split("-")[1]);

                        } else {
                            startRegion = Integer.parseInt(region.split("-")[0]);
                            endRegion = Integer.parseInt(region.split("-")[1]);
                        }
                    }
                    solenoidAnalysis.analysis(r.getPdbCode(), r.getPdbChain(), startRegion, endRegion);
                } catch (Exception ex) {

                    ex.printStackTrace();


                }

            }
        }
        // left-handed 1qre

        //
//        List<String> rows = DataIO.readLines("data/evaluation/solenoid.data.NoTRs");
//
//        for (String row : rows) {
//            String pdbCode = row.substring(0, 4);
//            String pdbChain = row.substring(5, 6);
//            solenoidAnalysis.analysis(pdbCode, pdbChain);
//        }

//        // right-handed
//        solenoidAnalysis.analysis("1bd8", "A");
//        solenoidAnalysis.analysis("1air", "A");
//        solenoidAnalysis.analysis("1pe9", "A");
//        // left-handed
//        solenoidAnalysis.analysis("1qre", "A");
//        solenoidAnalysis.analysis("1hm9", "A");
//        solenoidAnalysis.analysis("1j2z", "A");

//        solenoidAnalysis.demo();

    }

    int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {
            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;
        }
        return 0;
    }

    public void demo() throws StructureException {

//
        Vector3D o = new Vector3D(0, 0, 0);
        Vector3D x = new Vector3D(1, 0, 0);
        Vector3D y = new Vector3D(0, 1, 0);
        Vector3D t = new Vector3D(2, 1, 0);

        Plane plane = new Plane(o, x, y);
        System.out.println(plane.getNormal());
//
//
//        Vector3D vector1 = x.subtract(o);
//        Vector3D vector2 = t.subtract(o);
//        System.out.println(Math.toDegrees(Vector3D.angle(vector2, vector1)));

//        Vector2D o = new Vector2D(0, 0);
//        Vector2D x = new Vector2D(1, 0);
//        Vector2D y = new Vector2D(0, 1);
//        Vector2D t = new Vector2D(1, 1);
//        System.out.println(this.angleBetweenTwoPointsWithFixedPoint(t, y, o));

//        Atom[] atomS1 = new Atom[2];
//
//        atomS1[1] = new AtomImpl();
//        atomS1[1].setX(0);
//        atomS1[1].setY(0);
//        atomS1[1].setZ(0);
//        atomS1[0] = new AtomImpl();
//        atomS1[0].setX(1);
//        atomS1[0].setY(0);
//        atomS1[0].setZ(0);
//
//
//        Atom[] atomS2 = new Atom[2];
//        atomS2[1] = new AtomImpl();
//        atomS2[1].setX(0);
//        atomS2[1].setY(0);
//        atomS2[1].setZ(0);
//        atomS2[0] = new AtomImpl();
//        atomS2[0].setX(0);
//        atomS2[0].setY(1);
//        atomS2[0].setZ(0);
//
//        SVDSuperimposer svds = new SVDSuperimposer(atomS1, atomS2);
//        Matrix rotMatrix = svds.getRotation();
//        Atom tranMatrix = svds.getTranslation();
//        System.out.println(rotMatrix);
        //RotationAxis.getAngle(rotMatrix);

//        Superimposer superimpose = new Superimposer();
//        System.out.println(superimpose.superimposeSimple(atomS1,atomS2));
//        RotationAxis rotationAxis = new RotationAxis(rotMatrix, tranMatrix);
//        System.out.println(Math.toDegrees(rotationAxis.getAngle()));
//        boolean positiveScrew = Math.signum(rotationAxis.getRotationAxis().getX()) == Math.signum(rotationAxis.getScrewTranslation().getX());
//        System.out.println(positiveScrew);


    }


    public void analysis(String pdbCode, String pdbChain, int startRegion, int endRegion) throws StructureException {


        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);

        if (startRegion == 0 && endRegion == 0)
            endRegion = repeatFinder.getAtoms().length - 1;
        startRegion = this.getPosition(repeatFinder.getAtoms(), startRegion);
        endRegion = this.getPosition(repeatFinder.getAtoms(), endRegion);
        List<ProVector> lstVectorsOriignal = repeatFinder.getFeatures().getVectors();
        // extract six vectors
        lstVectorsOriignal = ProVector.toSecondaryVector(lstVectorsOriignal);
        List<ProVector> lstVectors = new ArrayList<ProVector>();
        for (ProVector vector : lstVectorsOriignal) {

            if (startRegion < vector.getPosStart() && vector.getPosStart() < endRegion)
                lstVectors.add(vector);

        }


        int sizeSegments = 6;
        int rotaion = 0;
        int size = 0;
        for (int i = 0; i < lstVectors.size() - sizeSegments; i++) {
            try {
                List<ProVector> subList = lstVectors.subList(i, i + sizeSegments);
                size++;
                rotaion += this.verifyRotation(subList);

            } catch (Exception ex) {


            }
            //break;
        }
        System.out.println(pdbCode + "_" + pdbChain + "\t" + rotaion);
    }

    public double verifyRotation(List<ProVector> lstVectors) throws StructureException {
        String output = "XXXXXXXXXXX";
        VectorModule vM = new VectorModule();

        List<Vector3D> points = new ArrayList<Vector3D>();
        for (ProVector vector : lstVectors) {
            points.add(new Vector3D(vector.getStart().getCoords()));
            points.add(new Vector3D(vector.getEnd().getCoords()));
        }


//        // find the axis rotation
//        int startX = 0;
//        int startY = 0;
//        int startZ = 0;
//        for (int i = 0; i < 3; i++) {
//            startX += points.get(i).getX();
//            startY += points.get(i).getY();
//            startZ += points.get(i).getZ();
//        }
//
//        Vector3D startLine = new Vector3D(startX / 3, startY / 3, startZ / 3);
//
//        startX = 0;
//        startY = 0;
//        startZ = 0;
//        for (int i = points.size() - 3; i < points.size(); i++) {
//            startX += points.get(i).getX();
//            startY += points.get(i).getY();
//            startZ += points.get(i).getZ();
//        }
//        Vector3D endLine = new Vector3D(startX / 3, startY / 3, startZ / 3);
//        Line line = new Line(startLine, endLine);
//
//
//        for (int i = 0; i < points.size() - 1; i++) {
//            Vector3D point1 = points.get(i);
//            Vector3D project1 = line.pointAt(line.getAbscissa(point1));
//
//            Vector3D point2 = points.get(i + 1);
//            Vector3D project2 = line.pointAt(line.getAbscissa(point2));
//
//            Atom[] atomS1 = new Atom[2];
//            atomS1[1] = new AtomImpl();
//            atomS1[1].setX(point1.getX());
//            atomS1[1].setY(point1.getY());
//            atomS1[1].setZ(point1.getZ());
//            atomS1[0] = new AtomImpl();
//            atomS1[0].setX(project1.getX());
//            atomS1[0].setY(project1.getY());
//            atomS1[0].setZ(project1.getZ());
//
//
//            Atom[] atomS2 = new Atom[2];
//            atomS2[1] = new AtomImpl();
//            atomS2[1].setX(point2.getX());
//            atomS2[1].setY(point2.getY());
//            atomS2[1].setZ(point2.getZ());
//            atomS2[0] = new AtomImpl();
//            atomS2[0].setX(project2.getX());
//            atomS2[0].setY(project2.getY());
//            atomS2[0].setZ(project2.getZ());
//
//            SVDSuperimposer svds = new SVDSuperimposer(atomS1, atomS2);
//            Matrix rotMatrix = svds.getRotation();
//            Atom tranMatrix = svds.getTranslation();
////            System.out.println(rotMatrix);
//            //RotationAxis.getAngle(rotMatrix);
//            RotationAxis rotationAxis = new RotationAxis(rotMatrix, tranMatrix);
//            System.out.println(Math.toDegrees(rotationAxis.getAngle()));
//            boolean positiveScrew = Math.signum(rotationAxis.getRotationAxis().getX()) == Math.signum(rotationAxis.getScrewTranslation().getX());
//            System.out.println(positiveScrew);
//            // find the point int the line
//            //Line line
//
//
//        }
//
//
//        for (Vector3D point : points) {
//            Vector3D projectPoint = line.pointAt(line.getAbscissa(point));
//
//        }


        Vector3D initialPoint = points.get(0);
        //find the next smallest point
        int indexNextPoint = 0;
        double distMin = 100;
        for (int i = 1; i < points.size(); i++) {
//            System.out.println(startPoint.distance(points.get(i)));
            double dist = initialPoint.distance(points.get(i));
            if (distMin > dist) {
                distMin = dist;
                indexNextPoint = i;
            }

        }

        // build a plane made by 3 points
        int numPointsInRound = indexNextPoint + 1;
        int L = (int) numPointsInRound / 3;
        if (L == 1)
            L = 1;
        Vector3D point1 = points.get(0);
        Vector3D point2 = points.get(L);
        Vector3D point3 = points.get(2 * L);
        Vector3D center3D = new Vector3D((point1.getX() + point2.getX() + point3.getX()) / 3, (point1.getY() + point2.getY() + point3.getY()) / 3, (point1.getZ() + point2.getZ() + point3.getZ()) / 3);


        Vector3D endPoint = points.get(points.size() - 1);


        Plane plane = new Plane(point1, point2, point3);
        Vector3D normalPlane = plane.getNormal();
        Vector3D direction = endPoint.subtract(center3D);
        //System.out.println(Math.toDegrees(Vector3D.angle(normalPlane,direction)));
        boolean isTheSamedirectionOfPlaneAndStructure = Math.toDegrees(Vector3D.angle(normalPlane, direction)) < 90;


        // projectPoint = plane.project(endPoint);


        // project point to plane;

        List<Vector2D> point2Ds = new ArrayList<Vector2D>();
//        for (int i = 0; i < numPointsInRound - 1; i++) {
//            Vector2D point2D = plane.toSubSpace(points.get(i));
//            point2Ds.add(point2D);
//        }

        point2Ds.add(plane.toSubSpace(point1));
        point2Ds.add(plane.toSubSpace(point2));
        point2Ds.add(plane.toSubSpace(point3));


//         calculate centre point

        double xCenter = 0.0;
        double yCenter = 0.0;
        for (Vector2D point : point2Ds) {
            xCenter += point.getX();
            yCenter += point.getY();

        }
        xCenter = xCenter / point2Ds.size();
        yCenter = yCenter / point2Ds.size();
        Vector2D center2D = new Vector2D(xCenter, yCenter);
        // choose one axis

        int rotation = 0;
        for (int i = 0; i < point2Ds.size() - 1; i++) {
            Vector2D startPoint = point2Ds.get(i);
            Vector2D nextPoint = point2Ds.get(i + 1);
//            System.out.println(this.angleBetweenTwoPointsWithFixedPoint(startPoint, nextPoint, center2D));
            double angle = this.angleBetweenTwoPointsWithFixedPoint(startPoint, nextPoint, center2D);
            if ((angle > 0 & isTheSamedirectionOfPlaneAndStructure) || (angle < 0 & !isTheSamedirectionOfPlaneAndStructure))
                rotation = rotation + 1;
            else
                rotation = rotation - 1;

        }


        //return this.angleBetweenTwoPointsWithFixedPoint(startPoint, nextPoint, center2D);

        return rotation;
    }

    Atom convertToAtom(Vector3D vector3D) {
        Atom atom = new AtomImpl();
        atom.setX(vector3D.getX());
        atom.setY(vector3D.getY());
        atom.setZ(vector3D.getZ());
        return atom;
    }

    private Vector3D project(Vector3D line1, Vector3D line2, Vector3D toProject) {
        double m = (double) (line2.getY() - line1.getY()) / (line2.getX() - line1.getX());
        double b = (double) line1.getY() - (m * line1.getX());

        double x = (m * toProject.getY() + toProject.getX() - m * b) / (m * m + 1);
        double y = (m * m * toProject.getY() + m * toProject.getX() + b) / (m * m + 1);

        return new Vector3D((int) x, (int) y);
    }

    public double angleBetweenTwoPointsWithFixedPoint(Vector2D point1,
                                                      Vector2D point2,
                                                      Vector2D fixedP) {

        double angle1 = Math.atan2(point1.getY() - fixedP.getY(), point1.getX() - fixedP.getX());
        double angle2 = Math.atan2(point2.getY() - fixedP.getY(), point2.getX() - fixedP.getX());
        double angle = Math.toDegrees(angle1 - angle2);
//        if (angle < 0) {
//            angle += 360;
//        }
        return angle;
    }

}
