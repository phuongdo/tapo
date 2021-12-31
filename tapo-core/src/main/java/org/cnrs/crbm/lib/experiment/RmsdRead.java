package org.cnrs.crbm.lib.experiment;

import org.cnrs.crbm.lib.io.DataIO;

import java.util.*;

/**
 * Created by pdoviet on 1/9/2015.
 */
public class RmsdRead {

    public static void main(String[] args) {


        Set<String> set = new HashSet<String>();
        List<String> lines = DataIO.readLines("data/rmsd.o");

        for (String line : lines) {

            if (line.length() > 0) {
                //System.out.println(line);
                set.add(line.substring(0, 6));
            }
        }

        System.out.println(set.size());
        for(String s: set){
            System.out.println(s);
        }
    }

}
