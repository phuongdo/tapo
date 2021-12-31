package org.cnrs.crbm.lib.io;

import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 4/17/2015.
 */
public class OutputReader {

    public static void main(String[] args) throws Exception {
        Map<String, List<Repeat>> output = new OutputReader().readResult("1lxa.o");
        System.out.println();

    }


    private Repeat convertToRepeat(String line) {
        Repeat repeat = new Repeat();
        String[] data = line.split("\t");
        try {
            String[] units = data[6].split(";");
            repeat.setScore(Double.parseDouble(data[7]));
            repeat.setRankScore(Double.parseDouble(data[8]));
            repeat.setFinderName(data[1]);
            for (String unit : units) {
                int start = 0;
                int end = 0;
                if (unit.startsWith("-")) {
                    unit = unit.substring(1, unit.length());
                    start = (-1) * Integer.parseInt(unit.split("-")[0]);
                    end = Integer.parseInt(unit.split("-")[1]);
                    repeat.getRepeats().add(new RepeatContent(start, end));

                } else {
                    start = Integer.parseInt(unit.split("-")[0]);
                    end = Integer.parseInt(unit.split("-")[1]);
                    repeat.getRepeats().add(new RepeatContent(start, end));
                }
            }

        } catch (Exception ex) {
        }

        return repeat;
    }


    public Map<String, List<Repeat>> readResult(String fileDir) throws Exception {

        //String fileDir = "data/tapo/ExperimentalResult.o";
        List<String> lines = DataIO.readLines(fileDir);
        Map<String, List<Repeat>> map = new HashMap<String, List<Repeat>>();
        for (String line : lines) {
            //System.out.println(line);
            if (line.length() < 10)
                continue;
            try {
                String[] data = line.split("\t");
                String pdbEntry = data[0];

                Repeat repeat = this.convertToRepeat(line);
                if (repeat.getScore() >= 0.0) {
                    if (map.containsKey(pdbEntry)) {
                        List<Repeat> l = map.get(pdbEntry);
                        l.add(repeat);
                        map.put(pdbEntry, l);
                    } else {
                        List<Repeat> l = new ArrayList<Repeat>();
                        l.add(repeat);
                        map.put(pdbEntry, l);
                    }
                }

            } catch (Exception ex) {
                //System.out.println();
            }
        }
        return map;

    }
}
