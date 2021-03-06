package org.cnrs.crbm.lib.otp.simulatedAnnealing;

import java.util.ArrayList;

public class TourManager {

	// Holds our cities
	private static ArrayList destinationCities = new ArrayList();

	// Adds a destination city
	public static void addCity(City city) {
		destinationCities.add(city);
	}

	// Get a city
	public static City getCity(int index) {
		return (City) destinationCities.get(index);
	}

	// Get the number of destination cities
	public static int numberOfCities() {
		return destinationCities.size();
	}

}