package org.cnrs.crbm.lib.repeats.module;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.geometry.euclidean.threed.Line;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.math.VectorCals;

public class VectorShape {

    public static void main(String[] args) {
        System.out.println(VectorShape.getSSPattern("BBBBB---HHHHHHHH--BBBB-HHH"));
    }

    /**
     * Get secondary pattern.
     *
     * @param secondaryStr
     * @return
     */
    public static String getSSPattern(String secondaryStr) {
        StringBuilder builder = new StringBuilder();
        int start = 0;
        int end = 0;
        for (int i = 0; i < secondaryStr.length(); i++) {
            char ch = secondaryStr.charAt(i);
            if (ch == 'B') {
                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'B') {
                    i++;
                }
                i = end = i - 1;
                // save
                if (end - start + 1 > 1)
                    builder.append('B');

            } else if (ch == 'H') {

                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'H') {
                    i++;
                }
                i = end = i - 1;

                if (end - start + 1 > 5)
                    builder.append('H');

            }

        }
        return builder.toString();
    }

    public static List<ProVector> getVectors(Atom[] atoms, String secondaryStr) {
        List<ProVector> vectors = new ArrayList<ProVector>();
        int start = 0;
        int end = 0;
        for (int i = 0; i < secondaryStr.length(); i++) {
            char ch = secondaryStr.charAt(i);
            if (ch == 'B') {
                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'B') {
                    i++;
                }
                i = end = i - 1;
                // save
                ProVector v = new ProVector();
                v.setStart(atoms[start]);
                v.setEnd(atoms[end]);
                v.setPosStart(start);
                v.setPosEnd(end);
                v.setOriginalStart(atoms[start]);
                v.setOriginalEnd(atoms[end]);
                v.setType("B");
                int size = end - start + 1;
                if (size > 2) {
                    List<ProVector> list = getProVector(atoms, start, end, 'B');
                    vectors.addAll(list);
                } else
                    vectors.add(v);

            } else if (ch == 'H') {

                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == 'H') {
                    i++;
                }
                i = end = i - 1;

                // check length first

                List<ProVector> list = getProVector(atoms, start, end, 'H');
                vectors.addAll(list);

            }

            // now , we consider to the loop regions!!!
            // updated from 14/04/2014
            else if (ch == '-') {
                start = i;
                while (i < secondaryStr.length()
                        && secondaryStr.charAt(i) == '-') {
                    i++;
                }
                i = end = i - 1;

                List<ProVector> list = getProVectorLoop(atoms, start, end);
                vectors.addAll(list);

            } // end loop regions

        }

        return vectors;
    }

    public static List<ProVector> getProVectorLoop(Atom[] atoms, int start,
                                                   int end) {
        // Sequence3D sequence3D = new Sequence3D();
        // String seq3d = "";
        // try {
        // String strAccs = "";
        // for (int i = 0; i < atoms.length; i++) {
        //
        // strAccs += "1";
        // }
        // seq3d = sequence3D.getSA(atoms, strAccs);
        //
        // } catch (Exception e) {
        // // TODO Auto-generated catch block
        // e.printStackTrace();
        // }

        // System.out.println(start + "-" + end + " : L");
        List<ProVector> vectors = new ArrayList<ProVector>();

        if (end - start + 1 >= 3) {
            for (int i = start; i <= end - 3; i = i + 3) {

                // System.out.println(seq3d.substring(i, i + 4));
                Vector3D v1 = new Vector3D(atoms[i].getCoords());
                Vector3D v2 = new Vector3D(atoms[i + 1].getCoords());
                Vector3D v3 = new Vector3D(atoms[i + 2].getCoords());
                Vector3D v4 = new Vector3D(atoms[i + 3].getCoords());
                Vector3D v12 = v2.subtract(v1);
                Vector3D v34 = v4.subtract(v3);
                ProVector v = new ProVector();
                double angle = Math.toDegrees(Vector3D.angle(v12, v34));
                if (angle > 120)
                    v.setType("m");
                else if (angle < 60)
                    v.setType("n");
                else
                    v.setType("l");

                Atom a_start = atoms[i];
                Atom a_end = atoms[i + 3];

                Atom clone_a_start = (Atom) a_start.clone();
                Atom clone_a_end = (Atom) a_end.clone();

                clone_a_start.setGroup(a_start.getGroup());
                clone_a_end.setGroup(a_end.getGroup());

                v.setPosStart(i);
                v.setPosEnd(i + 3);

                v.setStart(clone_a_start);
                v.setEnd(clone_a_end);
                v.setOriginalStart(a_start);
                v.setOriginalEnd(a_end);
                vectors.add(v);
            }
        } else {
            ProVector v = new ProVector();
            Atom a_start = atoms[start];
            Atom a_end = atoms[end];

            Atom clone_a_start = (Atom) a_start.clone();
            Atom clone_a_end = (Atom) a_end.clone();

            clone_a_start.setGroup(a_start.getGroup());
            clone_a_end.setGroup(a_end.getGroup());

            v.setPosStart(start);
            v.setPosEnd(end);
            v.setStart(clone_a_start);
            v.setEnd(clone_a_end);
            v.setType("x");
            vectors.add(v);

        }
        return vectors;

    }

    /**
     * convert secondary structure to a set of vectors helix : i and i + 4
     *
     * @param atoms
     * @param start
     * @param end
     * @return
     */
    public static List<ProVector> getProVector(Atom[] atoms, int start,
                                               int end, char type) {

        // System.out.println(start + "-" + end + " : " + type);

        // configure for type
        // H : i to i+4,
        // B : i to i+2

        // too long helix, dont break
        int LONG_HELIX = 100;

        int STEP = 4;
        double DIST_BREAK_VECT = 4.0;

        if (type == 'H') {
            STEP = 4;
            DIST_BREAK_VECT = 4.0;// no break;
        } else if (type == 'B') {
            STEP = 2;
            DIST_BREAK_VECT = 4.0;
        }

        List<ProVector> vectors = new ArrayList<ProVector>();

        Vector3D total_vector = new Vector3D(new double[]{0, 0, 0});

        for (int i = start; i <= end - STEP; i++) {

            Vector3D v = new Vector3D(atoms[i].getCoords());
            Vector3D v4 = new Vector3D(atoms[i + STEP].getCoords());
            // total_vector = total_vector + (v4 - v);
            total_vector = total_vector.add(v4.subtract(v));

        }

        Vector3D p1 = new Vector3D(atoms[start].getCoords());
        Vector3D p2 = new Vector3D(atoms[end].getCoords());
        p1 = VectorCals.project(p1, total_vector);
        p2 = VectorCals.project(p2, total_vector);
        double pa_distance = p1.distance(p2);
        Vector3D unit_v = total_vector.normalize();

        // calculate avg start point of vector
        if (STEP == 2)
            STEP = 3;
        p1 = new Vector3D(new double[]{0, 0, 0});
        for (int i = start; i < start + STEP - 1; i++) {
            Vector3D v = new Vector3D(atoms[i].getCoords());
            p1 = p1.add(v);

        }

        double a = (double) 1 / (STEP - 1);
        p1 = p1.scalarMultiply(a);
        p2 = p1.add(unit_v.scalarMultiply(pa_distance));

        // get the middle point?
        Line line = new Line(p1, p2);
        // middle point
        int middle = Math.round((end + start) / 2);

        Vector3D middle_point = new Vector3D(atoms[middle].getCoords());
        // distance from middle point to summary vector
        double m_dist = line.distance(middle_point);

        // System.out.println(m_dist);
        // /if > threshold, split this vector into 2 vectors

        Atom a_start = atoms[start];
        Atom a_end = atoms[end];
        Atom a_middle = atoms[middle];

        Atom clone_a_start = (Atom) a_start.clone();
        Atom clone_a_end = (Atom) a_end.clone();
        Atom clone_a_midle = (Atom) a_middle.clone();
        clone_a_start.setGroup(a_start.getGroup());
        clone_a_end.setGroup(a_end.getGroup());
        clone_a_midle.setGroup(a_middle.getGroup());
        if (m_dist > DIST_BREAK_VECT && ((end - start + 1) < LONG_HELIX)) {

            // break into 2 vectors here
            ProVector v1 = new ProVector();
            v1.setType(type + "");
            ProVector v2 = new ProVector();
            v2.setType(type + "");

            clone_a_start.setCoords(p1.toArray());
            clone_a_end.setCoords(p2.toArray());

            // set v1
            v1.setStart(clone_a_start);
            v1.setEnd(clone_a_midle);
            // set v2
            v2.setStart(clone_a_midle);
            v2.setEnd(clone_a_end);
            // set original if necessary

            v1.setOriginalStart(a_start);
            v1.setOriginalEnd(a_middle);

            v2.setOriginalStart(a_middle);
            v2.setOriginalEnd(a_end);

            v1.setPosStart(start);
            v1.setPosEnd(middle);
            v2.setPosStart(middle);
            v2.setPosEnd(end);

            vectors.add(v1);
            vectors.add(v2);

        } else {

            ProVector v = new ProVector();
            v.setType(type + "");

            clone_a_start.setCoords(p1.toArray());
            clone_a_end.setCoords(p2.toArray());

            v.setStart(clone_a_start);
            v.setEnd(clone_a_end);
            // set original if necessary
            v.setOriginalStart(a_start);
            v.setOriginalEnd(a_end);
            v.setPosStart(start);
            v.setPosEnd(end);
            vectors.add(v);

        }

        return vectors;
    }

}
