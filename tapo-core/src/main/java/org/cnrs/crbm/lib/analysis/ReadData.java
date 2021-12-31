package org.cnrs.crbm.lib.analysis;


import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.repeats.RepeatFinder;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.*;

/**
 * Created by pdoviet on 8/10/2015.
 */
public class ReadData {

    public static void main(String[] args) {

        List<String> lines = DataIO.readLines("data/francois/FrancoisRepeatSequence");
        HashMap<String, String> allFrancoisMap = new HashMap<String, String>();
        for (String line : lines) {
            String pdbCode = line.split("\t")[0];
            String repeatSeqs = line.split("\t")[1];
            allFrancoisMap.put(pdbCode, repeatSeqs);

        }

        lines = DataIO.readLines("data/benchmarkset/dataset3.in");
        Set<String> noFrancoisSet = new HashSet<String>();
        for (String line : lines) {
            String pdbCode = line.split("\t")[0];
            noFrancoisSet.add(pdbCode);

        }

        for (Map.Entry<String, String> entry : allFrancoisMap.entrySet()) {
            String key = entry.getKey();
            String value = entry.getValue();
            if (!noFrancoisSet.contains(key)) {
                System.out.println(key);
//                System.out.println(">" + key);
//                System.out.println(value);
//
//                System.out.println(">" + key);
//                String pdbCode = key.substring(0, 4);
//                String pdbChain = key.substring(5, 6);
//                RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
//                System.out.println(finder.getStrSeq());

            }


        }


    }


}
