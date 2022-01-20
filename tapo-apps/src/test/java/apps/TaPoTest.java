package apps;

import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.junit.jupiter.api.Test;

class TaPoTest {

    @Test
    void main() {
        RepeatFinder repeatFinder = new RepeatFinder("/input/1ib2.pdb", "1ib2", "A", 1);
        repeatFinder.findRepeat();
    }
}