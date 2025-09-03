package edu.put.ma.model.services;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import lombok.AllArgsConstructor;
import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@XmlRootElement(name = "strand")
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(propOrder = { "idx", "distance" })
@AllArgsConstructor
@NoArgsConstructor
@EqualsAndHashCode
public class DpStrand {
    @XmlElement
    @Getter
    @Setter
    private String idx;

    @XmlElement
    @Getter
    @Setter
    private String distance;

    @Override
    public String toString() {
        return new StringBuilder().append(idx).append(",").append(distance).toString();
    }
}
