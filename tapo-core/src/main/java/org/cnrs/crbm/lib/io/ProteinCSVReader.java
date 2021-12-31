package org.cnrs.crbm.lib.io;

import java.beans.PropertyDescriptor;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.bean.ColumnPositionMappingStrategy;
import au.com.bytecode.opencsv.bean.CsvToBean;

public class ProteinCSVReader {

	public List<Row> getData(String fileDir) throws Exception {
		CSVReader reader;
		reader = new CSVReader(new FileReader(fileDir), ',');
		String[] header = reader.readNext();
		CsvToBean<Row> bean = new CsvToBean<Row>();
		ColumnPositionMappingStrategy<Row> strategy = new ColumnPositionMappingStrategy<Row>();
		strategy.setType(Row.class);
		strategy.setColumnMapping(header);
		List<Row> rows = bean.parse(strategy, reader);
		// for (Row row : rows) {
		// System.out.println(row.getPdbCode() + ":" + row.getPdbChain());
		// }
		return rows;

	}

	public List<RowRepeatDB> getRepeatDB(String fileDir) throws Exception {
		CSVReader reader;
		reader = new CSVReader(new FileReader(fileDir), '\t');
		String[] header = reader.readNext();
		CsvToBean<RowRepeatDB> bean = new CsvToBean<RowRepeatDB>();
		ColumnPositionMappingStrategy<RowRepeatDB> strategy = new ColumnPositionMappingStrategy<RowRepeatDB>();
		strategy.setType(RowRepeatDB.class);
		strategy.setColumnMapping(header);
		List<RowRepeatDB> rows = bean.parse(strategy, reader);
		// for (Row row : rows) {
		// System.out.println(row.getPdbCode() + ":" + row.getPdbChain());
		// }

		return rows;

	}

}
