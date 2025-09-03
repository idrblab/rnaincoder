package edu.put.ma.model.services;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

import org.apache.commons.lang3.StringUtils;

import lombok.EqualsAndHashCode;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@XmlRootElement(name = "component")
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(propOrder = { "name", "description" })
@NoArgsConstructor
@EqualsAndHashCode
public class DrawComponent implements DpParameter {
    @XmlElement
    @Getter
    @Setter
    private String name;

    @XmlElement
    @Getter
    @Setter
    private String description;

    public DrawComponent(final String pname) {
        this.name = pname;
    }

    public DrawComponent(final String pname, final String pdescription) {
        this(pname);
        this.description = pdescription;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder().append('"').append(name);
        if (StringUtils.isNotBlank(description)) {
            sb.append(":").append(description);
        }
        return sb.append('"').toString();
    }
}
