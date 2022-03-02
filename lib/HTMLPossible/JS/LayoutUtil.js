// This file is created to contain functions that aid with
// creating an organized layout for a page.
// Functions explained in project's README.txt



function lUaddBasicLayout(DOM_object, basic_layout_d) {
/*
 *  This function takes an existing DOM object with a parent
 *      DOM object and defines its foundational layout
 *      pieces such as left, top, height, width
 * Args: 
 *      DOM_object is a Document Object Model object
 *      basic_layout_d is an object with the following keys:
 *          values_type: (str) "fixed" or "fractions"
 *          l: (num) left or null
 *          t: (num) top or null
 *          h: (num) height or null
 *          w: (num) width or null
 */

    let lft = null;
    let tp = null;
    let ht = null;
    let wd = null;
    if (basic_layout_d["values_type"] === "fractions") {
        checkFracValues(basic_layout_d);
        let parent_DOM_elem = DOM_object.parentElement;
        if (!(parent_DOM_elem === undefined ||  parent_DOM_elem === null)) {
            // Note: we subtract 1 from height and width because of a seeming extension that occurs.
            if (!(basic_layout_d["l"] === null)) {
                lft = parent_DOM_elem.clientWidth*basic_layout_d["l"] - 1;
            }
            if (!(basic_layout_d["t"] === null)) {
                tp = parent_DOM_elem.clientHeight*basic_layout_d["t"] - 1;
            }
            if (!(basic_layout_d["h"] === null)) {
                ht = parent_DOM_elem.clientHeight*basic_layout_d["h"]-1;
            }
            if (!(basic_layout_d["w"] === null)) {
                wd = parent_DOM_elem.clientWidth*basic_layout_d["w"]-1;
            }
        } else {
            throw "In add basic layout, the DOM_object has no parent element"
        }
    } else if (basic_layout_d["values_type"] === "fixed") {
        if (!(basic_layout_d["l"] === null)) {
            lft= basic_layout_d["l"];
        }
        if (!(basic_layout_d["t"] === null)) {
            tp = basic_layout_d["t"];
        }
        if (!(basic_layout_d["h"] === null)) {
            ht = basic_layout_d["h"];
        }
        if (!(basic_layout_d["w"] === null)) {
            wd = basic_layout_d["w"];
        }
    } else {
        throw "value type must be either 'fractions' or 'fixed'"
    }

    lUAddLayoutSizeLocToDOMObj(DOM_object, lft, tp, ht, wd);

    return DOM_object

}

function lUAddLayoutSizeLocToDOMObj(DOM_obj, lft, tp, ht, wd) {
    if (!(lft === null)) {
        DOM_obj.style.left = lft.toString() + "px";
    }
    if (!(tp === null)) {
        DOM_obj.style.top =  tp.toString() + "px";
    }
    if (!(ht === null)) {
        DOM_obj.style.height = ht.toString() + "px";
    }
    if (!(wd === null)) {
        DOM_obj.style.width = wd.toString() + "px";
    }
}



function lUAddStyleToDOMObj(DOM_obj, style_d) {
    // For each key in style_d object
    // we add that style to the DomObject
    style_vals = Object.keys(style_d);
    for (let i = 0; i < style_vals.length; i++) {
        key = style_vals[i];
        func_str = "DOM_obj.style." + key + ' = "' + style_d[key] + '";';
        eval(func_str); 
    }
}


function lUAddPropertiesToDOMObj(DOM_obj, property_d) {

        // For each key in property_d object
        //
        //     // we add that property to the DomObject
        //
        //         // e.g. property_d:
        //
        //             //   target: "_blank"
        
    let property_vals = Object.keys(property_d);
      
    for (let i = 0; i < property_vals.length; i++) {
      
        let key = property_vals[i];
    
        let func_str = "DOM_obj." + key + ' = "' + property_d[key] + '";';
    
        eval(func_str); 
      
    }
      
}


function lUAddElementToParent(DOM_obj, parent_id) {

    prnt_dobj = document.getElementById("#" + parent_id)
    prnt_dobj.appendChild(DOM_obj);
}

function lUCreateElementFromInfo(inp_obj) {
    /*
     *
     * Args:
     *     inp_obj:
     *         tag_type: (str)
     *         id_i: (obj)
     *              parent_id: (str)
     *              id: (str)
     *         size_loc_i: (obj) Contains location (height, etc. info)
     *         style_i: (obj) Contains stylistic info.
     *         [unq_prp]: (obj) Contains properties unique to this element.
     */


    id_i = inp_obj["id_i"];


    let parent_id = id_i["parent_id"];
    let parent_DOM_elem = null;
    if (parent_id != "body") {
        parent_DOM_elem = document.getElementById(id_i["parent_id"])
    } else {
        parent_DOM_elem = document.getElementsByTagName("BODY")[0];
    }

    if (parent_DOM_elem == null || parent_DOM_elem == undefined) {
        console.log("No parent Dom element")
        console.log(id_i);
        console.log(inp_obj);
        throw "No parent DOM object for object with id: " + id_i["id"];
    }

    let current_DOM_elem = null;
    if (inp_obj["tag_type"].toUpperCase() != "SVG") {
        current_DOM_elem = document.createElement(inp_obj["tag_type"]);
        current_DOM_elem.id = id_i["id"];


        parent_DOM_elem.appendChild(current_DOM_elem);
        lUaddBasicLayout(current_DOM_elem, inp_obj["size_loc_i"]);
        lUAddStyleToDOMObj(current_DOM_elem, inp_obj["style_i"]);
        if ("unq_prp" in inp_obj) {
            lUAddPropertiesToDOMObj(current_DOM_elem, inp_obj["unq_prp"])
        }
    } else {
        // this is an SVG element
        // height and width are 100%
        current_DOM_elem = create_d3_SVG_in_parent(
                        parent_id, id_i["id"], 
                        border=inp_obj["style_i"]["border"],
                        position=inp_obj["style_i"]["position"])
    }

    return current_DOM_elem;

}

function checkFracValues(basic_layout_d) {
 /* In this function we check if the values in the layout dict are indeed
  * fractions (between 0 and 1) 
 *      basic_layout_d is an object with the following keys:
 *          values_type: (str) "fixed" or "percentage"
 *          l: (f)
 *          t: (f)
 *          h: (f)
 *          w: (f)
 */
    let key_list = ["l","t","h","w"];
    key_list.forEach(function (item, index) {
        if (!(basic_layout_d[item] > -0.01 && basic_layout_d[item] < 1.01)){
            throw "Error with fraction value " + item + ": " + basic_layout_d[item].toString();
        }
    });
}


function prepInt(inp_i) {
 //    inp_i is an int or float with no decimal nonzero digits, value will be converted into
 //    string and have commas: 1000000 -> '1,000,000'
 //    OR inp_i is '?' or null or undefined
    
 //    converting floats into ints and ints into strings and strings into lists
 //

    let x = null;
    if (inp_i != '?' && inp_i != null && inp_i != undefined) {
        x = inp_i.toString().split("");
    } else {
        console.log("inp_i : " + inp_i.toString());
        return '?'
    }

    let op_str = ''
    while (x.length > 3) {
        let c_char = x[0];
        x = x.slice(1);
        op_str += c_char;
        if (x.length % 3 == 0) {
            op_str += ","
        }
    }

    op_str += x.join(""); 

    return op_str

}


function fracToPrcString(flt) {
    // flt is a fraction to be turned into percent and cut short to 3 decimals
    // returns string
    if (flt <= 1.01 && flt > -0.01) {
        flt = (flt*100).toString();
    } else {

        return "Percent Error? " + flt.toString()
    }

    // We split the float into before decimal and after
    let l = null
    if (flt.includes(".")) {
        l = flt.split(".");
    } else {
        l = [flt, "0"];
    }
    console.log(l)

    // We round the after-decimal part
    let ad = l[1].slice(0,2)
    if (parseInt(ad[-1]) > 4) {
        ad = (ad[0] + 1).toString();
    } else {
        ad = ad[0].toString();
    }

    let op = l[0] + "." + ad;

    return op
}


function create_d3_SVG_in_parent(parent_id, id, 
                        border = null, position = null,
                        width = "100%", height = "100%") {
    // All args are str    
    // We create an SVG within the object indicated by parent_id

    if (parent_id == undefined || id == undefined) {
        throw "Within create d3 SVG, we need id and parent_id"
    }

    let new_svg = d3.select("#" + parent_id).append("svg")
            .attr("id", id)
            .attr("width", width)
            .attr("height", height);

    if (border != null) {
            new_svg.attr("border", border);
    }
    if (position != null) {
            new_svg.attr("position", position);
    }

    new_svg.append("g");
    
    return new_svg
}


function createEntireShellFromShellDataObject(shell_data_obj) {
    // We make all objects inside Shell Data Object
    // shell_data_obj must contain key "lyt_vls"
    // which has only objects that follow the rules
    // of lyt_vls subwindows

    lyt_vls = shell_data_obj["lyt_vls"]

    if (lyt_vls == undefined) {
        throw_str = "To create an entire shell from shell data object" +
            " the shell data object must contain property lyt_vls.";
        throw throw_str
    } else {
        subwindow_keys = Object.keys(lyt_vls)
        for (let k=0; k < subwindow_keys.length; k++) {
            lUCreateElementFromInfo(lyt_vls[subwindow_keys[k]])
        }
    }


}

function setPageSize(scale_val) {
    //scale_val is a Number
    
    var scale = 'scale(' + scale_val.toString() + ')';
    document.body.style.webkitTransform =  scale;    // Chrome, Opera, Safari
    document.body.style.msTransform =   scale;       // IE 9
    document.body.style.transform = scale;     // General

}



function createTableWithRows(table_id, rows_info) {
    /*
     *
     * Args:
     *     table_id: str (DOM id of table element) 
     *     rows_info: list<row_info>
     *          row_info: list<cell_info, cell_info, ...>
     *              cell_info: Object
     *                  type: (str) ["link" or "text"]
     *                  inner_text: (str) The text inside the cell
     *                  id: (str) The DOM id
     *                  [color]: Color of text
     *                  [textDecoration]: How the text will look
     *                  
     *                  IF type is link:
     *                      [href]: Link to go to
     *                      [func]: The onclick function
     *                          [func_params]: If func, then func_params
     *                              must be present.
     *  
     *
     */

    // We get the table
    let tbl = document.getElementById(table_id);
    

    // We iterate through each scaffold, giving it an entry in the table
    for (let r = rows_info.length-1; r >= 0; r--) {

        // Create an empty <tr> element and add it to the 1st position of the table:
        let row = tbl.insertRow(0);
        row.style.border = "1px solid black";

        let c_row_info = rows_info[r];

        for (let i=0; i < c_row_info.length; i++) {

            // What does cell_info look like?
            let cell_info = c_row_info[i];
            // We insert this cell into the row
            let c_cell = row.insertCell(0);

            // 'type' can be 'link' or 'text' so far.
            if (cell_info["type"] == 'link') {
                new_link = document.createElement("a");
                new_link.innerHTML = cell_info['inner_text']

                if ('cursor' in cell_info) {
                    new_link.style.cursor = cell_info['cursor'];
                } else {
                    new_link.style.cursor = "pointer";
                }
                if ('textDecoration' in cell_info) {
                    new_link.style.textDecoration = cell_info['textDecoration'];
                } else {
                    new_link.style.textDecoration = "underline";
                }
                new_link.id = cell_info['id'];
                if ('color' in cell_info) {
                    new_link.style.color = cell_info['color'];
                }
                if ('func' in cell_info){
                    if (!('func_params' in cell_info)) {
                        throw_str = 'If func is key in cell_info, then func_params'
                        throw_str.concat(' must also be a key.')
                        throw throw_str;
                    }
                    new_link.onclick = function() {
                        cell_info['func'](cell_info['func_params'])
                    }
                }
                if ('href' in cell_info) {
                    new_link.href = cell_info['href'];
                }

                c_cell.appendChild(new_link);
            } else if (cell_info["type"] == 'text') {
                text_tag = document.createElement(cell_info['text_tag']);
                text_tag.innerHTML = cell_info['inner_text'];
                text_tag.id = cell_info['id'];
                if ('textDecoration' in cell_info) {
                    text_tag.style.textDecoration = cell_info['textDecoration'];
                } 
                if ('color' in cell_info) {
                    text_tag.style.color = cell_info['color'];
                }
                c_cell.appendChild(text_tag);
            } else {
                throw "cell info must be 'type' or 'text'"
            }
            //row.insertCell(i)
        }
    }
}


function printLoadingSign(div_id, create_var,
                          loading_text="Loading...",
                          font_size=25,  x=20, y=20,
                          loading_sign_id='loading-sign',
                          ) {
    /*
     * We add a Loading sign to the top of the svg
     *
     * Args:
     *  
     *      d3svg is already the svg object
     *      create_var: Number from [0,1,2] whether to add, replace, or remove the loading sign
     *      loading_text: (str) The text to be added to the load sign
     *      loading_sign_id: (str) DOM id of loading text
     */

    loading_div = document.getElementById(div_id)
    
    if (create_var == 0) {
        console.log("Loading - ")
        loading_div.innerHTML = loading_text
        //makeText(d3svg, 'bold', font_size, x, y, loading_text, 'black', loading_sign_id) 
    } else if (create_var == 1) { 
        if (loading_div != null) {
            // We increase Loading by a dot, or remove all dots
            let load_text = loading_div.innerHTML
            if (load_text.slice(-3) == '...') {
                loading_div.innerHTML = "Loading"
            } else {
                loading_div.innerHTML = load_text + ".";
            }
        }
    } else {
        if (loading_div != null) {
            console.log("Done Loading")
            loading_div.innerHTML = ""
        }
    }

}
