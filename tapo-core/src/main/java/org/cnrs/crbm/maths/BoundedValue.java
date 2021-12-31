package org.cnrs.crbm.maths;


import java.lang.Integer;
import java.util.Observable;

/**
* (c) DFKI GmbH 26.9.2002 
* <p>
* Orginal C source:
* Author Mike Jackson - University of Edinburgh - 1999-2001 
* <p>
* The <code>BoundedValue</code> class creates objects that manage an 
* integer variable whose value can change within given minimum and
* maximum bounds (initially <code>Integer.MIN_VALUE
* .. Integer.MAX_VALUE</code>). <code>BoundedValue</code> objects
* generate events when their values change
* @author Christoph Lauer 
* @version 1.0
*/
public class BoundedValue extends Observable {
    // Default minimum value that values in the BoundedValue can take
    private int _minimumValue = Integer.MIN_VALUE;
    // Default maximum value that values in the BoundedValue can take
    private int _maximumValue = Integer.MAX_VALUE;
    // Value
    private int _value = 0;

    /**
       Create a <code>BoundedValue</code> with the given default value
       <code>value</code> which can take a value equal to any
       allowable integer and generate an event with this value as
       argument.
    */
    public BoundedValue(int value) {
        setValue(value);
    }

    /**
       Create a <code>BoundedValue</code> with the given default value
       <code>value</code> which can take a value within the given
       bounds and generate an event with this value as argument.
    */
    public BoundedValue(int value, int minimumValue, int maximumValue) {
        setMinimum(minimumValue);
        setMaximum(maximumValue);
        setValue(value);
    }

    /**
    Set current value to be <code>value</code> - ensuring that
    bounds are respected - and generate an event with the new
    value as argument 
    */
    public void setValue(int value) {
        // Ensure bounds are respected
        _value =
            Math.max(Math.min(value, _maximumValue), _minimumValue);
        setChanged();
        notifyObservers(new Integer(value));
    }

    /** Get current value */
    public int getValue() {
        return _value;
    }

    /** Get minimum value */
    public int getMinimumValue() {
        return _minimumValue;
    }

    /** Get maximum value */
    public int getMaximumValue() {
        return _maximumValue;
    }

    /** Get range of values the <code>BoundedValue</code> can take */
    public int getRange() {
        return _maximumValue - _minimumValue;
    }

    /**
    Decrement value - respecting minimum bound - and generate an 
    event with the new value as argument.
    */
    public void decrementValue() {
        setValue(--_value);
    }

    /**
    Increment value - respecting maximum bound - and generate an 
    event with the new value as argument.
    */
    public void incrementValue() {
        setValue(++_value);
    }

    /**
    Set minimum bound of allowable values to be
    <code>minimumValue</code> - possibly generating an event as
    the current value may need to be adjusted to respect the new
    bounds
    */
    public void setMinimum(int minimumValue) {
        _minimumValue = minimumValue;
        ensureCorrectBounds();
    }

    /**
    Set maximum bound of allowable values to be
    <code>maximumValue</code> - possibly generating an event as
    the current value may need to be adjusted to respect the new
    bounds
    */
    public void setMaximum(int maximumValue) {
        _maximumValue = maximumValue;
        ensureCorrectBounds();
    }

    // Ensure maximum bound >= minimum bound
    private void ensureCorrectBounds() {
        // Ensure _maximumValue >= _minimumValue
        if (_maximumValue < _minimumValue) {
            int temp = _minimumValue;
            _minimumValue = _maximumValue;
            _maximumValue = temp;
        }
        // Update to ensure _value respects new bounds
        setValue(_value);
    }
}
