package org.cnrs.crbm.lib.trsfinder;

/**
 * Created by pdoviet on 8/19/2015.
 */
public class Region {

    protected int start;
    protected int end;

    public Region(int start, int end) {
        super();
        this.start = start;
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int size() {
        return this.end - this.start + 1;
    }

    @Override
    public String toString() {

        return this.start + "-" + this.end;
    }
}
