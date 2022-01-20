package demo;

import org.cnrs.crbm.lib.repeats.RepeatFinder;

public class TapoMainTest {
    public static void main(String[] args) {
        RepeatFinder repeatFinder = new RepeatFinder("input/1ib2.pdb", "1ib2", "A", 1);
        repeatFinder.findRepeat();
    }
}
