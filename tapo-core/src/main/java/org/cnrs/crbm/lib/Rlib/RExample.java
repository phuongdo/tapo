package org.cnrs.crbm.lib.Rlib;

/**
 * Created by pdoviet on 10/3/2014.
 */


import java.awt.Container;
import javax.swing.*;
import rcaller.RCaller;

public class RExample extends Container {


    ImageIcon image;


    public RExample() {

        RCaller caller = new RCaller();

        // Path of your R localisation
        caller.setRScriptExecutableFile("F:/App/R-3.1.0/bin/Rscript.exe");
        StringBuffer code = new StringBuffer();

        int i = 100;
        code.append("x<-rnorm(" + i + ", 0, 2);");
        code.append("png('C:/Users/pdoviet/Desktop/file.png');");
        code.append("plot.ts(x, xlab=\"pos\", ylab=\"contact\");");
        try {
            caller.RunRCode(code.toString(), false, false);
        } catch (Exception e) {
            e.printStackTrace();
        }
        // code.append("install.packages('adk'),");
        // code.append("library(spatstat);");
        // code.append("mydata <- read.csv('C:/Users/nznassi/Documents/test.csv', header = TRUE, sep=';');");
        // code.append("x <- mydata$x;");
        // code.append("y <- mydata$y;");
        // code.append("X <- ppp(x, y);");
        // code.append("png('C:/Users/nznassi/Pictures/file.png');");
        // code.append("plot(X);");
        // code.append("plot(density.ppp(X));");



    }

    public static void main(String[] args) {
        new RExample();
    }



}
