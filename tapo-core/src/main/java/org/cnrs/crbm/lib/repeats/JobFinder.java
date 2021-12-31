package org.cnrs.crbm.lib.repeats;

import org.cnrs.crbm.lib.trsfinder.Finder;

import java.util.concurrent.Callable;

/**
 * Created by pdoviet on 11/21/2014.
 */
public class JobFinder implements Callable {
    private Finder finder = null;

    public JobFinder(Finder finder) {
        this.finder = finder;
    }


    public Object call() {
        // starting process
        finder.start();
        return finder;
    }
}
