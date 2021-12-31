package org.cnrs.crbm.lib.utils;

import java.awt.*;
import java.awt.event.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public final class ColorUtils {

    private static ColorUtils instance = null;
    Map<Integer, Color> colors = new HashMap<Integer, Color>();

    private ColorUtils() {
        // /Random randomGenerator = new Random();

        // colors.put(0, Color.yellow);
        // colors.put(1, Color.green);
        // colors.put(2, Color.red);
        // colors.put(3, Color.blue);
        // colors.put(4, Color.orange);
        // colors.put(5, Color.black);
        // colors.put(6, Color.magenta);
        // colors.put(7, Color.cyan);
        // colors.put(8, Color.yellow);
//		colors.put(0, new Color(192, 208, 255));
//		colors.put(1, new Color(176, 255, 176));
//		colors.put(2, new Color(255, 192, 200));
//		colors.put(3, new Color(255, 255, 128));
//		colors.put(4, new Color(255, 192, 255));
//		colors.put(5, new Color(176, 240, 240));
//		colors.put(6, new Color(255, 208, 112));
//		colors.put(7, new Color(240, 128, 128));
//		colors.put(8, new Color(245, 222, 179));
//		colors.put(9, new Color(0, 191, 255));
//		colors.put(10, new Color(205, 92, 92));
//		colors.put(11, new Color(102, 205, 170));
//		colors.put(12, new Color(154, 205, 50));
//		colors.put(13, new Color(238, 130, 238));
//		colors.put(14, new Color(0, 206, 209));
//		colors.put(15, new Color(0, 255, 127));
//		colors.put(16, new Color(60, 179, 113));
//		colors.put(17, new Color(0, 0, 139));
//		colors.put(18, new Color(189, 183, 107));
//		colors.put(19, new Color(0, 100, 0));
//		colors.put(20, new Color(128, 0, 0));
//		colors.put(21, new Color(128, 128, 0));
//		colors.put(22, new Color(128, 0, 128));
//		colors.put(23, new Color(0, 128, 128));
//		colors.put(24, new Color(184, 134, 11));
//		colors.put(25, new Color(128, 0, 0));
//		colors.put(26, new Color(178, 34, 34));

        colors.put(0, new Color(60, 153, 245));
        colors.put(1, new Color(234, 65, 65));

    }

    public final static ColorUtils getInstance() {
        // Le "Double-Checked Singleton"/"Singleton doublement vérifié" permet
        // d'éviter un appel coûteux à synchronized,
        // une fois que l'instanciation est faite.
        if (ColorUtils.instance == null) {
            // Le mot-clé synchronized sur ce bloc empêche toute instanciation
            // multiple même par différents "threads".
            // Il est TRES important.

            ColorUtils.instance = new ColorUtils();

        }
        return ColorUtils.instance;
    }

    public Map<Integer, Color> getColors() {
        return this.colors;
    }

}