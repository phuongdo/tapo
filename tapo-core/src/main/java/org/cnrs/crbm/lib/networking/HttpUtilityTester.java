package org.cnrs.crbm.lib.networking;

/**
 * Created by pdoviet on 2/5/2015.
 */

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class HttpUtilityTester {

    /**
     * This program uses the HttpUtility class to send a GET request to
     * Google home page; and send a POST request to Gmail login page.
     */
    public static void main(String[] args) {


        // test sending POST request
        Map<String, String> params = new HashMap<String, String>();
        String requestURL = "http://protein.bio.unipd.it/raphael/Raphael.jsp";
        params.put("pdbcode", "1lxa");
        params.put("chain", "A");

        try {
            HttpUtility.sendPostRequest(requestURL, params);
            String[] response = HttpUtility.readMultipleLinesRespone();
            for (String line : response) {
                System.out.println(line);
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        HttpUtility.disconnect();
    }
}