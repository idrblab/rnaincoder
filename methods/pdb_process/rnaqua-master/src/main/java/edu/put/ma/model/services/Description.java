package edu.put.ma.model.services;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import lombok.Getter;
import lombok.Setter;

@XmlRootElement(name = "basics")
@XmlType(propOrder = { "filename", "errorList" })
public class Description {
    @Getter
    @Setter
    private String filename = null;

    @XmlElementWrapper(name = "errors")
    @XmlElement(name = "error")
    private List<String> errorList = new ArrayList<String>();

    public final List<String> getErrorLists() {
        return errorList;
    }

    public final void setErrorList(final List<String> perrorList) {
        errorList = perrorList;
    }

    public final void addAllErrorList(final List<String> perrorList) {
        errorList.addAll(perrorList);
    }

    public final void addError(final String error) {
        errorList.add(error);
    }

    public final int getErrorsCount() {
        return errorList.size();
    }

    @Override
    public final String toString() {
        final StringBuilder res = new StringBuilder(filename).append("\n");
        for (int i = 0; i < errorList.size(); i++) {
            res.append(errorList.get(i)).append("\n");
        }
        return res.toString();
    }
}
