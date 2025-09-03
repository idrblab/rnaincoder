package edu.put.ma;

import edu.put.ma.CommandEnum.Mapping;

public interface Command {

    Mapping getMapping();

    String toString();
}
