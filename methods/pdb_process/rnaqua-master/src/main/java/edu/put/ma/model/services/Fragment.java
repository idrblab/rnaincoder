package edu.put.ma.model.services;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import lombok.Getter;
import lombok.Setter;

@XmlRootElement(name = "fragment")
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(propOrder = { "residue", "distance", "resNo" })
public class Fragment {
    @XmlElement
    @Getter
    @Setter
    protected String residue;

    @XmlElement
    @Getter
    @Setter
    protected String distance;

    @XmlElement
    @Getter
    @Setter
    protected String resNo;
}
