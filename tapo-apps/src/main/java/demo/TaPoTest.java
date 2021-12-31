package demo;

import org.cnrs.crbm.lib.repeats.FinderMode;
import org.cnrs.crbm.lib.repeats.RepeatFinder;

/**
 * Created by pdoviet on 4/17/2015.
 */
public class TaPoTest {
    public static void main(String[] args) {
        // Get current time
        long start = System.currentTimeMillis();

        // Do something ...
        // 2yo3, 1mr7
        String pdbCode = "2yo3";
        String pdbChain = "A";
        int nCore = 6;

//        RepeatFinder repeatFinder = new RepeatFinder("input/rosmann_model.pdb", "pdbId", " ", nCore);
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain, nCore);
//        System.out.println(repeatFinder.getStrSS());
        repeatFinder.setMode(FinderMode.CONSOLE);
        System.out.println("starting...");
        repeatFinder.findRepeat();
        System.out.println(repeatFinder.getPdbCode() + "_" + repeatFinder.getPdbChain());
        System.out.println("PDB Length:" + repeatFinder.getAtoms().length);
        System.out.println("Repeat: " + repeatFinder.isRepeat());
        System.out.println(repeatFinder.getTAPOScore());
        System.out.println(repeatFinder.getConsoleOutputDetails());

        // Get elapsed time in milliseconds
        long elapsedTimeMillis = System.currentTimeMillis() - start;

        // Get elapsed time in seconds
        float elapsedTimeSec = elapsedTimeMillis / 1000F;


        // Get elapsed time in minutes
        float elapsedTimeMin = elapsedTimeMillis / (60 * 1000F);

        // Get elapsed time in hours
        float elapsedTimeHour = elapsedTimeMillis / (60 * 60 * 1000F);

        // Get elapsed time in days
        float elapsedTimeDay = elapsedTimeMillis / (24 * 60 * 60 * 1000F);

        System.out.println(elapsedTimeSec);

    }
}
