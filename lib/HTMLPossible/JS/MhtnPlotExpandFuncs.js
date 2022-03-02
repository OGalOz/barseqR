
/*
 *
 *    Process: First we create all the DOM objects, then we populate them.
 *
         
    INFO:

    Functional Programming Approach. Data is in file ./TestMhtnDefaults.js
        name of data is 'MhtnDefaults'.
    //SECTIONS: What about Return to All Scaffolds Button?
    // FullDisplayDiv
    //    GRAPH DISPLAY DIV
    //         Graph Title Div
    //         Graph SVG
    //         Graph_Info_Sidebar
    //              Point_Info_DIV
    //                  Point Info Title
    //                  InfoScrollDownText
    //              Scaffold_List_DIV
    //                  Scaffolds_List_Title
    //                  Scaffold_List_ScrollDown
    //              All_Scaffolds_View_Button
    //    EXPLANATION DIV
    //         Explanation Text
   
        
    // FUNCTIONS
    // CreateMhtnDisplay
    //     CreateGraphDisplay
    //         CreateGraphAndTitleDiv
    //              CreateGraphTitleDiv
    //              CreateGraphSVGDiv
    //                  CreateGraphSVG
    //                      PopulateGraphSVG
    //         CreateInfoSideBar
    //              CreatePointInfoDiv
    //                  CreatePointInfoTitle
    //                  CreateInfoScrollDownText
    //              CreateScaffoldListDiv
    //                  CreateScaffoldsListTitle
    //                  CreateScaffoldListScrollDown
    //              CreateAllScaffoldsViewButton
    //     CreateExplainDisplayDiv
    //         CreateExplainText
    //
    //         
    //      PopulateGraphSVG
    //          MakeAllScaffoldsGraph
    //              MakeAllScaffoldAxes
    //                  MakeLine
    //                  MakeXAxisLabel
    //                      MakeText
    //


 
    Input Data:
        Mhtn_data:
            'scaffolds' -> scaffolds_data_obj (d)
                scaffolds_data_obj:
                    scaffold_name -> scaffold_obj
                        scaffold_obj: 
                            scaffold_length: (int)
                            max_z: (float)
                            pos_to_Zscr_l: list<[pos, Z_scr]>
            'mean': float
            'SD': float
            'max_z': float
            'analysis_type' : (str) "AllGenomeStats" or "IndividualScaffoldStats" 
             

            

*/


function CreateMhtnDisplay(MhtnDefaults, MH_start_data, MH_data) {
        /*
         *
         * Args:
         *     MhtnDefaults:
         *         object that contains the defaults for the Manhattan plot
         *     MH_start_data:
         *         object that contains variables for the Manhattan plot, e.g.
         *            which div contains the visualization?
         *     MH_data:
         *         object that contains the data for the visualization
         *
         */

    // We maintain an object relating object Tag IDs to their DOM values.
    // Note this already exists in document.getElementById function
    const id_to_DOM_obj_d = {};

    CreateGraphDisplay(MhtnDefaults, MH_start_data, MH_data);
    CreateExplainDisplayDiv(MhtnDefaults, MH_start_data);
    PopulateElements(MH_data, MhtnDefaults, MH_start_data);
}

function PopulateElements(MH_data, MhtnDefaults, MH_start_data) {
    /*
     * In this segment of the code, we take existing elements and populate them with content.
     * Upon loading, we populate using the "All Scaffolds" Values
     */

    let svg_d3 = d3.select("#" + MhtnDefaults["lyt_vls"]["graph_svg"]["id_i"]["id"]);
    let scaffolds_info_and_maxSD = OrderScaffoldNamesAndLengths(MH_data, svg_d3); 
    // scfld_length_array: list<[scaffold_name (str), scaffold_length (int), scaffold_start (int)]>
    let scf_names_and_lengths = scaffolds_info_and_maxSD[0]
    // max_z_for_all_scf: Number: the total max Z-Score
    let max_z = scaffolds_info_and_maxSD[1]
    // Do we want the following map to be global?
    let scf_to_color_obj = CreateScaffoldToColorMap(scf_names_and_lengths, MhtnDefaults)

    PopulateScaffoldInfoTable(MhtnDefaults, scf_names_and_lengths, 
                              scf_to_color_obj, svg_d3, max_z);

    PopulateGraphSVG(MhtnDefaults, MH_start_data, MH_data, 
                     scf_names_and_lengths, max_z, scf_to_color_obj)



}

function CreateScaffoldToColorMap(scf_names_and_lengths, MhtnDefaults) {

    let scf_to_color_obj = {};
    let color_list = MhtnDefaults["Scaffold_Info_Table_Values"]["color_options"];

    for (let i = 0; i < scf_names_and_lengths.length; i++) {
        scf_to_color_obj[scf_names_and_lengths[i][0]] = color_list[i % color_list.length]
    }

    return scf_to_color_obj
}

function PopulateScaffoldInfoTable(MhtnDefaults, scf_names_and_lengths, scf_to_color_obj, svg_d3, max_z) {
    /*
     *
     * Args:
     *      MhtnDefaults: The giant Default data object passed in
     *      scf_names_and_lengths: list<list< scaffold_name (str), scaffold_length (int), scaffold_start (int) >>
     *      scf_to_color_obj: map from scaffold name to color name
     *      svg_d3: The d3 svg object
     *
     *      
     *
     *
     */

    // We get the table
    let tbl = document.getElementById(MhtnDefaults["lyt_vls"]["scfld_list_tbl"]["id_i"]["id"]);
    
    // We shorten the name of the variable
    let sc_l = scf_names_and_lengths;


    // We iterate through each scaffold, giving it an entry in the table
    for (let r = sc_l.length-1; r >= 0; r--) {

            // Create an empty <tr> element and add it to the 1st position of the table:
            let row = tbl.insertRow(0);
            row.style.border = "1px solid black";
            
            // Insert new cells (<td> elements) at the 1st and 2nd position of the "new" 
            // <tr> element:
            let cell1 = row.insertCell(0);
            let cell2 = row.insertCell(1);
            let cell3 = row.insertCell(2);
            
            // Add scaffold name and then scaffold length
            
            // First the ranking of scaffold in terms of length
            cell1.innerHTML = (r + 1).toString();
            // Then the scaffold name and link
            let scf_name = sc_l[r][0];
            console.log(scf_name);
            let scf_link = document.createElement("a");
            scf_link.innerHTML = scf_name;
            scf_link.style.cursor = "pointer";
            scf_link.style.textDecoration = "underline";
            scf_link.style.color = scf_to_color_obj[scf_name];
            // Here is where we set the on click function to open up that single scaffold's graph
            scf_link.onclick = function() {
                MakeSingleScaffoldGraph(scf_name, MH_data, scf_to_color_obj, sc_l, max_z )
            }
            scf_link.id = ScaffoldNameToValidHTMLID(scf_name) + "-link";
            cell2.appendChild(scf_link);
            cell3.innerHTML = prep_int(sc_l[r][1]) + " bp";
        }

}



function CreateGraphDisplay(MhtnDefaults, MH_start_data, MH_data) {
    let graph_display_elem = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["graph_display_div_i"]);
    CreateGraphAndTitleDiv(MhtnDefaults, MH_start_data)
    CreateInfoSideBar(MhtnDefaults, MH_start_data, MH_data)
}

function CreateGraphAndTitleDiv(MhtnDefaults, MH_start_data) {

    let graph_and_title_div_elem = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["graph_and_title_div"]);
    CreateGraphTitleDiv(MhtnDefaults, MH_start_data)
    CreateGraphSVGDiv(MhtnDefaults, MH_start_data, MH_data)

}

function CreateGraphTitleDiv(MhtnDefaults, MH_start_data) {

    let graph_title_elem = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["graph_title_div"]);
    graph_title_elem.innerHTML = MhtnDefaults["lyt_vls"]["graph_title_div"]["graph_title"]
}

function CreateGraphSVGDiv(MhtnDefaults, MH_start_data, MH_data) {

    let graph_svg_div_elem = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["graph_svg_div"]);
    CreateGraphSVG(MhtnDefaults, MH_start_data, MH_data);
}

function CreateGraphSVG(MhtnDefaults, MH_start_data, MH_data) {
    svg_info = MhtnDefaults["lyt_vls"]["graph_svg"];

    let svg = d3.select("#" + svg_info["id_i"]["parent_id"]).append("svg")
            .attr("width", "100%")
            .attr("height", "100%" )
            .attr("border", svg_info["style_i"]["border"])
            .attr("position", svg_info["style_i"]["position"])
            .attr("id", svg_info["id_i"]["id"])
            .append("g");
    
    return svg

}

function CreateInfoSideBar(MhtnDefaults, MH_start_data, MH_data) {
    let info_side_bar = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["graph_info_sidebar"]);
    CreatePointInfoDiv(MhtnDefaults, MH_start_data, MH_data);
    CreateScaffoldListDiv(MhtnDefaults, MH_start_data, MH_data);
    CreateAllScaffoldsViewButton(MhtnDefaults, MH_start_data, MH_data);
}

function CreatePointInfoDiv(MhtnDefaults, MH_start_data, MH_data) {
    let point_info_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["point_info_div"]);
    CreatePointInfoTitle(MhtnDefaults, MH_start_data, MH_data);
    CreateInfoScrollDownText(MhtnDefaults, MH_start_data, MH_data);
}

function CreatePointInfoTitle(MhtnDefaults, MH_start_data, MH_data) {

    let point_info_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["point_info_ttl_div"]);
    point_info_div.innerHTML = "Point Info:"

}

function CreateInfoScrollDownText(MhtnDefaults, MH_start_data, MH_data) {
    
    let point_info_scrldwn_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["point_info_scrldwn_div"]);

}

function CreateScaffoldListDiv(MhtnDefaults, MH_start_data, MH_data) {

    let scfld_list_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["scfld_list_div"]);
    CreateScaffoldsListTitle(MhtnDefaults, MH_start_data, MH_data);
    CreateScaffoldListScrollDown(MhtnDefaults, MH_start_data, MH_data);
}

function CreateScaffoldsListTitle(MhtnDefaults, MH_start_data, MH_data) {
    let scfld_list_ttl_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["scfld_list_ttl_div"]);
    scfld_list_ttl_div.innerHTML = "Scaffolds:"
}

function CreateScaffoldListScrollDown(MhtnDefaults, MH_start_data, MH_data) {

    let scfld_list_scrldwn_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["scfld_list_scrldwn_div"]);

    // The below function creates a table in some location 
    CreateScaffoldInfoTable(MhtnDefaults) 
}

function CreateAllScaffoldsViewButton(MhtnDefaults, MH_start_data, MH_data) {

    let scfld_list_scrldwn_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["all_scaffolds_btn_div"]);

}

function CreateExplainDisplayDiv(MhtnDefaults, MH_start_data) {
    let explain_display_div = LUCreateElementFromInfo(
            MhtnDefaults["lyt_vls"]["explanation_div"]);
    explain_display_div.innerHTML = MhtnDefaults["lyt_vls"]["explanation_div"]["explain_text"].join("");
}


function PopulateGraphSVG(MhtnDefaults, MH_start_data, MH_data, scf_names_and_lengths, max_z, scf_to_color_obj) {
    // This is the most complex function in the program
    MakeAllScaffoldsGraph(MH_data, MhtnDefaults, scf_names_and_lengths, max_z, scf_to_color_obj);
}


function GetDefaultVariables() {
        x_i = MhtnDefaults["x_i"];
        y_i = MhtnDefaults["y_i"];
}
           

function SetPageSize(scale_val) {
    //scale_val is a Number
    
    var scale = 'scale(' + scale_val.toString() + ')';
    document.body.style.webkitTransform =  scale;    // Chrome, Opera, Safari
    document.body.style.msTransform =   scale;       // IE 9
    document.body.style.transform = scale;     // General

}

function RemoveReturnButton() {

        //Remove return button if exists
        let return_btn = document.getElementById("return-btn");
        if (return_btn) {
            return_btn.parentNode.removeChild(return_btn);
        }
}

function MakeAllScaffoldsGraph(MH_data, MhtnDefaults, scf_names_and_lengths, max_z, scf_to_color_obj) {
    /*
     * MH_data: The entire data for the scaffolds etc.
     * MhtnDefaults: imported from TestMhtnDefaults.js
     *
        // scfld_length_array: list<[scaffold_name (str), scaffold_length (int), scaffold_start (int)]>
     *
        // max_z_for_all_scf: Number: the total max Z-Score
     */

    // Change d3 select to ID
    //

        // Clear old SVG
        let svg_obj = document.getElementById(MhtnDefaults["lyt_vls"]["graph_svg"]["id_i"]["id"]);
        let svg_d3 = d3.select("#" + svg_obj.id);
        svg_d3.selectAll("*").remove();

        // If the return button exists, we remove it.
        RemoveReturnButton();


        let y_ticks = GetProperTicks(0, Math.round(max_z + 1));
        let x_ticks = scf_names_and_lengths;
        let x_i = MhtnDefaults["Mhtn_Display_SVG_Values"]["x_i"];
        let y_i = MhtnDefaults["Mhtn_Display_SVG_Values"]["y_i"];
        let org = MhtnDefaults["Mhtn_Display_SVG_Values"]["org"];
        let svg_i = {
            "width": svg_obj.clientWidth,
            "height": svg_obj.clientHeight,
            "scaffold_bool": true
        }

        let z = MakeAllScaffoldAxes( svg_d3, org, x_ticks, y_ticks, svg_i, x_i, y_i);

        let x_axis_len = z[0];
        let y_axis_len = z[1];
        let total_bp_len = z[2];
        let x_org = z[3];
        let y_org = z[4];
        // Create these from data (1% is max, .2% is min?)
        let max_crc_r = svg_i["width"]*MhtnDefaults[
                                    "Mhtn_Display_SVG_Values"]["GraphValues"]["max_circle_frac"];
        let min_crc_r = svg_i["width"]*MhtnDefaults[
                                    "Mhtn_Display_SVG_Values"]["GraphValues"]["min_circle_frac"];


        let graph_info = {
            "org": [x_org, y_org], 
            "svg_d3": svg_d3,
            "min_circle_radius": min_crc_r,
            "max_circle_radius": max_crc_r,
            "x_axis_len": x_axis_len,
            "y_axis_len": y_axis_len,
            "max_y_val": y_ticks[y_ticks.length - 1],
            "total_bp_len": total_bp_len,
            "scf_to_color": scf_to_color_obj 
        };


        PopulateManhattanSkyline(scf_names_and_lengths, 
            MH_data["scaffolds"], graph_info,
            false, MH_data); 
        

}


function AddPositionDiv() {
        // This function adds the div and text which describes the point info
       
        p_i = MhtnDefaults["Position_Div_Info"]

        h_i = p_i["h_i"]
        // Adding the Info div
        var h = document.createElement(h_i["tag_type"]); // Create the H1 element 
        h.style.left = h_i["left"];
        h.style.top = h_i["top"];
        h.innerHTML = h_i["base_text"];
        h.id = h_i["h_id"]


        var pos_info_dv = document.createElement("div");
        pos_info_dv.style.top = p_i["top"];
        pos_info_dv.style.width = p_i["width"];
        pos_info_dv.style.position = p_i["position"];
        pos_info_dv.id = p_i["div_id"];
        pos_info_dv.appendChild(h);
        document.body.appendChild(pos_info_dv)

}

//HERE
function AddReturnButton(scf_names_and_lengths, max_z, scf_to_color_obj) {

        let rtn_btn = document.getElementById("return-btn");

        if (!rtn_btn) {

            let rb_d = MhtnDefaults["lyt_vls"]["all_scf_btn"];

            let return_btn_obj = LUCreateElementFromInfo(rb_d);
            return_btn_obj.innerHTML = rb_d["innerHTML"];
            return_btn_obj.onmouseover = function () {
                return_btn_obj.style.backgroundColor = rb_d["hover_color"]; 
            }
            return_btn_obj.onmouseout = function () {

                return_btn_obj.style.backgroundColor = rb_d["style_i"]["backgroundColor"] 
            }
            return_btn_obj.onclick = function() {MakeAllScaffoldsGraph(
                    MH_data, MhtnDefaults, scf_names_and_lengths, max_z, scf_to_color_obj
                    )};
            /*
            var return_btn = document.createElement("button");
            return_btn.style.top = rb["top"];
            return_btn.style.left = rb["left"];
            return_btn.style.width = rb["width"];
            return_btn.style.height = rb["height"];
            return_btn.style.position = rb["position"] ;
            return_btn.style.paddingBottom = rb["paddingBottom"];
            return_btn.style.backgroundColor = rb["bg_color"];
            return_btn.innerHTML = '<h3 style="color:' + rb["text_color"] + '">' + rb["inner_text"] + '</h3>';
            return_btn.id = rb["btn_id"];
            return_btn.style.cursor = "pointer";
            return_btn.onclick = function() {MakeAllScaffoldsGraph()};
            document.body.appendChild(return_btn)
            */

        }


}


function AddTitle(org, x_axis_len) {

        // Args: org: [x_origin, y_origin], x_axis_len: (Number)

        // Adding the title
        var h = document.createElement("H1"); // Create the H1 element 
        h.style.left = (org[0] + x_axis_len/2 - 30).toString() + "px";
        h.style.top = "0px";
        h.style.position = "absolute";
        h.innerHTML = "Insertions Plot";
        // var t = document.createTextNode(); // Create a text element 
        // h.appendChild(t); // Append the text node to the H1 element 
        document.body.appendChild(h); // Append the H1 element to the document body 
}




function PopulateManhattanSkyline(scf_names_and_lengths, scfs_info, graph_info,
                                    single_scaff_bool, MH_data) {
    /*
     * Args:
     *  scf_names_and_lengths: list<[scaffold_name (str), scaffold_length (int), scaffold_start (int)]>
     *  scfs_info:
     *      scaffold_name -> scaffold_info_d
     *          scaffold_info_d: (d) Contains keys
     *              scaffold_length: (str)
     *              max_z: (Number)
     *              mean: (Number)
     *              SD: (Number)
     *              pos_to_Zscr_l: list<[position, SD above mean]>
     *              [color]: string (If included retains the color of the scaffold)
     *  graph_info: (d)
     *      org: [x Number, y Number] Graph origin
     *      svg_d3: d3 object
     *      min_circle_radius: (Number)
     *      max_circle_radius: (Number)
     *      x_axis_len: Number
     *      y_axis_len: Number
     *      max_y_val: (Number) Represents highest Z-Score ,
     *          All circle sizes are filled in ratio to this number
     *      scf_to_color: Object mapping scaffold name to color
     *  single_scaff_bool: (Boolean) Whether or not it's displaying a single scaffold
     *     total_bp_len: Number (sum of lengths of all scaffolds)
     */

    // Here we maintain a list of scaffolds to colors for later expansion
    let next_scf_x_start_point = graph_info["org"][0];

    for ( i=0; i < scf_names_and_lengths.length; i++) {
        let scf_name = scf_names_and_lengths[i][0]
        let scaffold_info_d = scfs_info[scf_name];

        // We do a deep copy of the object
        new_scf_i_d = JSON.parse(JSON.stringify(scaffold_info_d));
        
        graph_info["x_scaffold_start"] = next_scf_x_start_point;
        graph_info["current_color"] = graph_info["scf_to_color"][scf_name]; 


        if (MH_data["analysis_type"] == "AllGenomeStats") {
            new_scf_i_d["mean"] = MH_data["mean"];
            new_scf_i_d["SD"] = MH_data["SD"];
            new_scf_i_d["max_z"] = MH_data["max_z"];
        }

        next_scf_x_start_point = CreateCirclesForScaffold(new_scf_i_d, 
                                graph_info, scf_names_and_lengths[i][0]);

    }


}


function CreateCirclesForScaffold(new_scf_i_d, graph_info, scaffold_name) {
    /*
     *
     * new_scf_i_d: (d) Contains keys
     *     scaffold_length: (int)
     *     mean: (Number)
     *     SD: (Number)
     *     max_z: (Number)
     *     pos_to_Zscr_l: list<[position, SD above mean]>
     *
     * graph_info: (d)
     *     svg_d3: d3 object
     *     min_circle_radius: (Number)
     *     max_circle_radius: (Number)
     *     max_y_val: (Number) Represents highest z score. 
     *          All circle sizes are filled in ratio to this number
     *     current_color: (string) The color of these circles in this scaffold
     *     org: [x Number, y Number] Graph origin
     *     x_axis_len: Number
     *     y_axis_len: Number
     *     x_scaffold_start: Number
     *     total_bp_len: Number (sum of lengths of all scaffolds)
     *
     *
     * scaffold_name: (str) 
     */

    let pos_l = new_scf_i_d["pos_to_Zscr_l"];

    //console.log(pos_l)

    let minc = graph_info["min_circle_radius"];
    let maxc = graph_info["max_circle_radius"];
   
    // BELOW NEEDS TO BE UPDATED 
    let scf_x_len = (new_scf_i_d["scaffold_length"]/graph_info["total_bp_len"])*(
                                                            graph_info["x_axis_len"]);

    for (let j=0; j <pos_l.length; j++ ) {
        // pos_SD is a list of [position (int), SD (number)]
        let pos_SD = pos_l[j];

        // Calculating the center of the circle. We only draw circles with SD > 0;
        if (pos_SD[1] > 0 && new_scf_i_d["scaffold_length"] > 0 && graph_info["max_y_val"] > 0) {
            let circle_radius =  minc + (maxc - minc)*(pos_SD[1]/graph_info["max_y_val"]);
            let cx = graph_info["x_scaffold_start"] + (pos_SD[0]/new_scf_i_d["scaffold_length"])*scf_x_len
            let cy = graph_info["org"][1] - graph_info["y_axis_len"]*(pos_SD[1]/graph_info["max_y_val"])

            // Drawing the circle
            graph_info["svg_d3"].append("circle")
                                .attr("cx", cx)
                                .attr("cy", cy)
                                .attr("r", circle_radius)
                                .attr("fill", graph_info["current_color"])
                                .attr("opacity", 0.5)
                                .data([pos_SD])
                                .on("click", function(d) {
                                    UpdatePointInfo(d, scaffold_name, new_scf_i_d)
                                    
                                });

        }

    }

    // Returning the next starting point
    return graph_info["x_scaffold_start"] + scf_x_len;

}

function UpdatePointInfo(point_info_array, scaffold_name, new_scf_i_d) {
    /* Args
     *      point_info_array: [position, SD above mean]
     *
     */

    pos_obj = document.getElementById(MhtnDefaults[
            "lyt_vls"]["point_info_scrldwn_div"]["id_i"]["id"]);
    inner_HTML = "<p>Scaffold: " + scaffold_name + "</p> "; 
    point_pos = point_info_array[0].toString();
    inner_HTML += "<p>Position: " + prep_int(point_pos) + "</p>";
    num_insertions = Math.round(
        new_scf_i_d["mean"] + new_scf_i_d["SD"]*point_info_array[1]
    );
    inner_HTML += "<p> # Insertions: " + prep_int(num_insertions) + "</p>";
    pos_obj.innerHTML = inner_HTML;

}

function GetSortedScaffoldNamesAndLengths(MH_data) {
    /*
     * This function sorts all the scaffolds in decreasing length.
     *
     *  inp_d: 
     *      analysis_type: (str)
     *      scaffolds: scaffold_info_d
     *          scaffold_info_d:
     *              scaffold_name (str) -> scaffold_info_d (d)
     *                  scaffold_info_d: (d) Contains keys
     *                     scaffold_length: (int)
     *
     *
     * Returns
     *
     *      sorted_scf_name_length_array: list<list<scf_name, scf_length>>
     */

    let scfld_info_d = MH_data["scaffolds"]
    let scfld_names = Object.keys(scfld_info_d);

    // We sort scaffolds by length of scaffold
    let sorted_scf_name_length_array = []
    for (let i = 0; i<scfld_names.length; i++) {
        let c_scf = scfld_names[i]; 
        sorted_scf_name_length_array.push([c_scf, scfld_info_d[c_scf]["scaffold_length"]])
    }
    sorted_scf_name_length_array.sort(SortByNum);

    return sorted_scf_name_length_array

}

function OrderScaffoldNamesAndLengths(inp_d, svg_obj) {
    /*
     * This function sorts all the scaffolds in decreasing length.
     *
     *  inp_d: 
     *      analysis_type: (str)
     *      scaffolds: scaffold_info_d
     *          scaffold_info_d:
     *              scaffold_name (str) -> scaffold_info_d (d)
     *                  scaffold_info_d: (d) Contains keys
     *                     scaffold_length: (int)
     *                     max_z: (Number)
     *                     mean: (Number)
     *                     SD: (Number)
     *                     pos_to_Zscr_l: list<[position, SD above mean]>
     *
     * Returns
     *      tuple: <scfld_length_array, max_z_for_all_scf>
     *          scfld_length_array: list<[scaffold_name (str), scaffold_length (int), scaffold_start (int)]>
     *          max_z_for_all_scf: Number: the total max Z-Score 
     *
     */


    sorted_scf_name_length_array = GetSortedScaffoldNamesAndLengths(inp_d)
    
    let scfld_names = [];
    for (let i=0; i<sorted_scf_name_length_array.length; i++) {
        scfld_names.push(sorted_scf_name_length_array[i][0])
    }


    
    let scfld_info_d = inp_d["scaffolds"];
    let scfld_length_array = [];
    let max_z_for_all_scf = 0;
    let total_bp_scf_start = 0;

    for (let i = 0; i<scfld_names.length; i++) {

        scfld_length_array.push([scfld_names[i], scfld_info_d[scfld_names[i]]["scaffold_length"],
                                total_bp_scf_start]);
        total_bp_scf_start += scfld_info_d[scfld_names[i]]["scaffold_length"];

        if (inp_d["analysis_type"] == "IndividualScaffoldStats") {
            if (scfld_info_d[scfld_names[i]]["max_z"] > max_z_for_all_scf) {
                max_z_for_all_scf = scfld_info_d[scfld_names[i]]["max_z"];
            }
        }
    };

    if (inp_d["analysis_type"] == "AllGenomeStats") {
        max_z_for_all_scf = inp_d["max_z"]
    }

    return [scfld_length_array, max_z_for_all_scf];

}

function CreateScaffoldInfoTable(MhtnDefaults) {
    /* We create the scaffold info table inside the div 'scfld-list-scrldwn-div'
     *
     */


    /* The point of the function is to make a table in javascript
     * which lists the scaffolds in order and their lengths.
     * 
     *
     * Args:
     *     sorted_scf_names_and_lengths: list<scf_tp>
     *         scf_tp: [scf_length (int), scf_name (str)]
     * 
     */
    
    // We get sorted scaffold names and lengths

    let scf_tbl_obj = LUCreateElementFromInfo(MhtnDefaults["lyt_vls"]["scfld_list_tbl"]);
    let scf_tbl_entries_i = MhtnDefaults["Scaffold_Info_Table_Values"]["entry_size_info"];
}


function ScaffoldNameToValidHTMLID(scf_name) {

    newscf = scf_name.replace("/","-");

    return newscf

}


function SortByNum(a,b) {
    if (a[1] < b[1]) return 1;
   if (a[1] > b[1]) return -1;
   return 0;
}


function MakeAndGetSVG(parent_div_id, width, height, left, right, top, bottom) {
    /*
     * Args:
     *      parent_div_id: (str)
     *      width -> the rest: Number describing Svg
     */

        var svg = d3.select("#" + parent_div_id).append("svg")
            .attr("width", width)
            .attr("height", height )
            .attr("border", "1px solid #F0F0F0")
            .attr("transform", 
                  "translate(" + left + "," + top + ")")
            //.call(d3.zoom().on("zoom", function() {
            //    svg.attr("transform", d3.event.transform)
            //})
            //)
            .append("g");
        return svg

}

function MakeLine(svg_obj, color, x1, y1, x2, y2, stroke_width ) {
    /*
     * Args: 
     *  svg_obj: A d3 svg object
     *  color: str, like "black"
     *  x1 - y2, Numbers
     *  stroke_width: width of line, Number
     *
     * Note: We need to make sure the numbers are relatively close to integers
     */

               return svg_obj.append('line')
                   .attr('x1', x1.toFixed(2))
                   .attr('y1', y1.toFixed(2))
                   .attr('x2', x2.toFixed(2))
                   .attr('y2', y2.toFixed(2))
                   .attr('stroke', color)
                   .attr('stroke-width', stroke_width)
                   .attr('position', 'absolute')
                   ;

}

function MakeText(svg_obj, font_weight, font_size, x, y, text_str) {
    /*
     *  Args:
     *  
     *      svg_obj: A d3 svg object
     *      font_weight: (str) like "bold", "normal",
     *      font_size: Number
     *      x, y: Number
     *      text_str: (str) Text you want to make
     *
     */
            svg_obj.append('text')
                .attr('font-weight', font_weight)
                .attr('font-size', font_size)
                .attr('x', x)
                .attr('y', y)
                .text(text_str);

}

//function MakeRect() {};

//function MakeCircle() {};
function MakeAllScaffoldAxes(svg_d3, org, x_ticks, y_ticks, svg_i, x_i, y_i) {
            /*
             *
             * TO DO: Make creation of axes function based.
             * Make plus button for each scaffold
            Here we create the axes - 
            Args:
               svg_d3: the d3 svg object to which we append values
               org: Graph origin (tuple) <x_start (FRACTION), y_start (FRACTION)>
               x_ticks: list<[scaffold_name (str), scaffold_length (int), scaffold_start (int)]>
               y_ticks: list<int> The numbers in the y axis
               x_i: x_axis info - comes from MhtnDefaults
                    x_title_i:
                        size_loc_i
                        label
                        style_i
                    x_ticks_i:
                        size_loc_i
                        tick_length: (Number) fraction of enclosing box which is tick length
                        style_i
                    x_axis_i:
                        size_loc_i
                        style_i
               y_i: y_axis_info - comes from MhtnDefaults
                    y_title_i:
                        size_loc_i
                        label
                        style_i
                    y_ticks_i:
                        size_loc_i
                        tick_length: (Number) fraction of enclosing box which is tick length
                        style_i
                    y_axis_i:
                        size_loc_i
                        style_i
                svg_i:
                    width: Number
                    height: Number
                    scaffold_bool: (Boolean)
            Returns:
                [x_axis_len (Number), y_axis_len (Number), single_x_tick_len (Number)]
            */
            // The following (x_ticks) is a list with tuples of length 3 (IF FULL GENOME)
            // It is a list of ints if a single scaffold


            let Xmin_num = 0;
            // We calculate lengths of all scaffolds together:
            Xmax_num = 0;
            for (i=0; i< x_ticks.length; i++) {
                Xmax_num += x_ticks[i][1]
            }

            let Ymin_num = y_ticks[0];
            let Ymax_num = y_ticks[y_ticks.length -1];
            let Ydist = Ymax_num - Ymin_num;
            
            // We create the axis lengths 
            x_axis_len = x_i["x_axis_i"]["size_loc_i"]["w"]*(svg_i["width"]);
            y_axis_len = y_i["y_axis_i"]["size_loc_i"]["h"]*(svg_i["height"]);

            // x origin
            // y origin
            x_org = org[0]*svg_i["width"];
            y_org = org[1]*svg_i["height"];

            
            //Making X axis line (no ticks)
            MakeLine(svg_d3, "black", x_org, y_org, x_org + x_axis_len, y_org, 
                 x_i["x_axis_i"]["style_i"]["strokeWidth"]);

            //Making Y axis line (no ticks)
            MakeLine(svg_d3, "black", x_org, y_org, x_org, y_org - y_axis_len, 
                 y_i["y_axis_i"]["style_i"]["strokeWidth"]);
            

            // Labels
            // X-Axis Text Label
            MakeXAxisLabel(svg_d3, svg_i, x_i["x_title_i"], x_i["x_title_i"]["label"]);
            //MakeText(svg_d3, "bold", x_i["x_label_font_size"], org[0] + x_axis_len/3,
            //        org[1] + x_i["x_label_dst"] + x_i["x_tick_len"], x_i["x_label"])
           
            MakeYAxisLabel(svg_d3, svg_i, y_i["y_title_i"], y_i["y_title_i"]["label"]);

            Add_X_Ticks_To_SVG_All_Scf(svg_d3, x_org, y_org, x_i, Xmax_num, x_axis_len, x_ticks, svg_i);
            Add_Y_Ticks_To_SVG(svg_d3, x_org, y_org, y_i, Ydist, y_axis_len, Ymin_num, y_ticks, svg_i);

            return [x_axis_len, y_axis_len, Xmax_num, x_org, y_org]

}

function MakeXAxisLabel(svg_d3, svg_i, title_i, title_str) {
     // We place the text about 1/3 down the X Axis
     // 2/3 Down the location of Y
   

    let x = svg_i["width"] * (title_i["size_loc_i"]["l"] + title_i["size_loc_i"]["w"] * (1/3));
    let y = svg_i["height"] * (title_i["size_loc_i"]["t"] +  title_i["size_loc_i"]["h"] * (2/3))

    MakeText(svg_d3, title_i["style_i"]["fontWeight"],
            title_i["style_i"]["fontSize"],
            x, y, 
            title_str)

}
function MakeYAxisLabel(svg_d3, svg_i, title_i, title_str) {


            Yx = svg_i["width"] * (title_i["size_loc_i"]["l"] + title_i["size_loc_i"]["w"] * (1/2))
            Yy = svg_i["height"] * (title_i["size_loc_i"]["t"] +  title_i["size_loc_i"]["h"] * (1/2))

            ax_tsl = Yx.toString() + "," + Yy.toString()

            svg_d3.append('text')
                .attr('font-weight', title_i["style_i"]["fontWeight"])
                .attr("transform", "translate(" + ax_tsl + ") rotate(270)")
                .attr('font-size', title_i["style_i"]["fontSize"])
                .text(title_str);

}

function Add_X_Ticks_To_SVG_Single_Scf(svg_d3, x_org, y_org, x_i, Xmax_num, x_axis_len, x_ticks, svg_i, Kb_bool) {
        /*
         * x_ticks: list<num>
         * y_ticks: list<num>
         *
         */

            let tks_i = x_i["x_ticks_i"]
            let xtik_len = svg_i["height"]*tks_i["size_loc_i"]["h"]*tks_i["tick_length"];
            let dist_to_txt = svg_i["height"]*tks_i["size_loc_i"]["h"]*tks_i["distance_to_text"];
   
            for (let i=0; i < x_ticks.length; i++) {
                let crnt_tick_info = x_ticks[i];
                
                let x_loc = x_org + (crnt_tick_info/Xmax_num)*(x_axis_len)

                // Add the x tick  
                MakeLine(svg_d3, tks_i["style_i"]["color"], x_loc, y_org, x_loc, 
                            y_org + xtik_len, tks_i["style_i"]["strokeWidth"]);
               
                if (x_ticks.length > 1 && i > 0) {
                        if (!(i == x_ticks.length - 2)) {
                            //Here we make it so the second to last tick has no text
                            // for spacing purposes.
                            if (Kb_bool) {
                               c_tick_label_s = prep_int(Math.round(crnt_tick_info/1000));
                            } else {
                                c_tick_label_s = prep_int(crnt_tick_info); 
                            }

                            MakeText(svg_d3, tks_i["style_i"]["fontWeight"], tks_i["style_i"]["fontSize"], 
                     x_loc - tks_i["width_shift"], y_org + xtik_len + dist_to_txt,
                            c_tick_label_s);
                            
                            //MakeText(svg_d3, "normal", x_i["x_ticks_font_size"], 
                             //     x_loc, org[1] + x_i["x_tick_len"] + 15,
                             //     c_tick_label_s);
                        }

                }


            }


            //Add the last tick for x-axis
            let x_loc = x_org + x_axis_len;
            MakeLine(svg_d3, tks_i["style_i"]["color"], x_loc, y_org, x_loc, 
                            y_org + xtik_len, tks_i["style_i"]["strokeWidth"]);

            /*MakeLine(svg_d3, tks_i["style_i"]["color"], x_loc, y_org, x_loc, 
                            y_org + tick_len, tks_i["style_i"]["strokeWidth"]);*/

}


function Add_X_Ticks_To_SVG_All_Scf(svg_d3, x_org, y_org, x_i, Xmax_num, x_axis_len, x_ticks, svg_i) {

            let x_loc = x_org 
            let tks_i = x_i["x_ticks_i"]
    
            let tick_len = svg_i["height"]*tks_i["size_loc_i"]["h"]*tks_i["tick_length"]; 

            for (let i=0; i < x_ticks.length; i++) {
                let crnt_tick_info = x_ticks[i];

                // Add the x tick  
                MakeLine(svg_d3, tks_i["style_i"]["color"], x_loc, y_org, x_loc, 
                            y_org + tick_len, tks_i["style_i"]["strokeWidth"]);

                x_loc += (crnt_tick_info[1]/Xmax_num)*(x_axis_len);
                
            }
            //Add the last tick for x-axis
            MakeLine(svg_d3, tks_i["style_i"]["color"], x_loc, y_org, x_loc, 
                            y_org + tick_len, tks_i["style_i"]["strokeWidth"]);

}

function Add_Y_Ticks_To_SVG(svg_d3, x_org, y_org, y_i, Ydist, y_axis_len, Ymin_num, y_ticks, svg_i) {

            let tik_i = y_i["y_ticks_i"];
            let ytik_len = svg_i["width"]*tik_i["size_loc_i"]["w"]*tik_i["tick_length"];
            let dist_to_txt = svg_i["width"]*tik_i["size_loc_i"]["w"]*tik_i["distance_to_text"];

            for (let i=0; i < y_ticks.length; i++) {
                ytick = y_ticks[i];
                // We get y location
                let y_loc = y_org - y_axis_len*(ytick - Ymin_num)/(Ydist);

                //Add the y tick  
                MakeLine(svg_d3, tik_i["style_i"]["color"], x_org, y_loc, x_org - ytik_len,
                            y_loc,  tik_i["style_i"]["strokeWidth"]);

                MakeText(svg_d3, tik_i["style_i"]["fontWeight"], tik_i["style_i"]["fontSize"], 
                        x_org - ytik_len - dist_to_txt, y_loc + tik_i["height_shift"],
                        ytick.toString());
            }
}

function CreateToolTip() {



    //We create a tooltip to describe the y axis

    tool_text = "Explanation: We calculate the mean by taking the total number of insertions ^" +
                "  and dividing that by the total number of locations with insertions.^" +
                "We compute the Standard Deviation using that mean. The Z-score for a ^" +
                "point is (#Insertions - Mean)/SD. We do not list any values whose Z-score ^" +
                "is below 0.";


    var b = tool_text.split('^');

    for (i = 0; i < b.length; i++) {
        c_txt_s = b[i];
        text_elem = document.createElement("P");
        text_elem.innerHTML = c_txt_s;
        text_elem.style.position = "absolute";
        text_elem.style.left = "400px";
        text_elem.style.top = (760 + i*20).toString() + "px";
        text_elem.style.fontWeight = "900"
        document.body.appendChild(text_elem)
    }

}


function MakeSingleScaffoldAxes(svg_d3, org, x_ticks, y_ticks, svg_i, x_i, y_i) {
            /*
             *
             * TO DO: Make creation of axes function based.
             * Make plus button for each scaffold
             * NOTE: There is a difference between circle location and tick location.
             *  Ticks are spread out evenly, whereas circles aren't. More calculations
             *  necessary.
            Here we create the axes - 
            Args:
               svg_d3: the d3 svg object to which we append values
               org: Graph origin (tuple) <x_start (Number), y_start (Number)>
               x_ticks: list<num>
               y_ticks: list<num>
               x_i: x_axis info - comes from MhtnDefaults
                    x_title_i:
                        size_loc_i
                        label
                        style_i
                    x_ticks_i:
                        size_loc_i
                        tick_length: (Number) fraction of enclosing box which is tick length
                        style_i
                    x_axis_i:
                        size_loc_i
                        style_i
               y_i: y_axis_info - comes from MhtnDefaults
                    y_title_i:
                        size_loc_i
                        label
                        style_i
                    y_ticks_i:
                        size_loc_i
                        tick_length: (Number) fraction of enclosing box which is tick length
                        style_i
                    y_axis_i:
                        size_loc_i
                        style_i
               svg_i:
               width: Number
               height: Number
               scaffold_bool: (Boolean)
            Returns:
                [x_axis_len (Number), y_axis_len (Number), single_x_tick_len (Number)]
            */
            // The following (x_ticks) is a list with tuples of length 3 (IF FULL GENOME)
            // It is a list of ints if a single scaffold
            let Xmin_num = 0;
            // The final X value is the start of the last scaffold plus its length
            let Xmax_num = x_ticks[x_ticks.length - 1]

            let Kb_bool = false;
            if (x_ticks.length > 1) {

                if (x_ticks[1] - x_ticks[0] > 1000) {
                    Kb_bool = true;
                }

            }

            let Ymin_num = y_ticks[0];
            let Ymax_num = y_ticks[y_ticks.length -1];
            let Ydist = Ymax_num - Ymin_num;
            
            // We create the axis lengths 
            let x_axis_len = x_i["x_axis_i"]["size_loc_i"]["w"]*(svg_i["width"]);
            let y_axis_len = y_i["y_axis_i"]["size_loc_i"]["h"]*(svg_i["height"]);

           
            // x origin
            // y origin
            let x_org = org[0]*svg_i["width"];
            let y_org = org[1]*svg_i["height"];



            //Making X axis line (no ticks)
            MakeLine(svg_d3, "black", x_org, y_org, x_org + x_axis_len, y_org, 
                 x_i["x_axis_i"]["style_i"]["strokeWidth"]);

            //Making Y axis line (no ticks)
            MakeLine(svg_d3, "black", x_org, y_org, x_org, y_org - y_axis_len, 
                 y_i["y_axis_i"]["style_i"]["strokeWidth"]);
            

            let x_axis_label_s = "?" 
            if (Kb_bool) {
                x_axis_label_s = x_i["x_title_i"]["label"] + " Location (Kbp)"
            } else {
                x_axis_label_s = x_i["x_title_i"]["label"] + " Location (bp)"
            }
            
            // Labels
            // X-Axis Text Label
            MakeXAxisLabel(svg_d3, svg_i, x_i["x_title_i"],x_axis_label_s);
            MakeYAxisLabel(svg_d3, svg_i, y_i["y_title_i"], y_i["y_title_i"]["label"]);

          
            Add_X_Ticks_To_SVG_Single_Scf(svg_d3, x_org, y_org, x_i, Xmax_num, x_axis_len, x_ticks, svg_i, Kb_bool);

            


            Add_Y_Ticks_To_SVG(svg_d3, x_org, y_org, y_i, Ydist, y_axis_len, Ymin_num, y_ticks, svg_i);


            return [x_axis_len, y_axis_len, Xmax_num]

}  



//NOW
function MakeSingleScaffoldGraph(scf_name, MH_data, scf_to_color_obj, scf_names_and_lengths, max_z) {
    /* This function embodies the entire process of creating a graph
     * Args:
     *      svg_obj: d3 SVG object
     *      scf_name: (str) Name of scaffold
     *
     */

    // We redraw the entire graph - clear SVG
    // Redo Y and X axis with single scaffold
    // X axis will have no ticks
    // Need 'go back' button.
    // Change X axis label to be just Scaffold Name

    // Clear old SVG

    let svg_obj = document.getElementById(MhtnDefaults["lyt_vls"]["graph_svg"]["id_i"]["id"]);
    let svg_d3 = d3.select("#" + svg_obj.id);
    svg_d3.selectAll("*").remove();

    /*
    d3.selectAll("svg > *").remove();
    svg_obj = d3.select("svg")
    */

    //Make return button
    AddReturnButton( scf_names_and_lengths, max_z, scf_to_color_obj);

    // scf_info_d contains keys: scaffold_length, pos_to_Zscr_l
    let current_scf_info_d = MH_data["scaffolds"][scf_name]
    let current_y_max_z =  Math.round(current_scf_info_d["max_z"])
    let x_ticks = GetProperTicks(0, current_scf_info_d["scaffold_length"]);
    let y_ticks = GetProperTicks(0, Math.round(current_scf_info_d["max_z"] + 1))

    let x_i = MhtnDefaults["Mhtn_Display_SVG_Values"]["x_i"];
    let y_i = MhtnDefaults["Mhtn_Display_SVG_Values"]["y_i"];
    let org = MhtnDefaults["Mhtn_Display_SVG_Values"]["org"];
    let svg_i = {
        "width": svg_obj.clientWidth,
        "height": svg_obj.clientHeight,
        "scaffold_bool": true
    }

    let z = MakeSingleScaffoldAxes(svg_d3, org, x_ticks, y_ticks, svg_i, x_i, y_i);

    //NOW WE DRAW THE ACTUAL CIRCLES

    let x_axis_len = z[0];
    let y_axis_len = z[1];
    let X_max_num = z[2];

    let max_crc_r = svg_i["width"]*MhtnDefaults[
                                "Mhtn_Display_SVG_Values"]["GraphValues"]["max_circle_frac"];
    let min_crc_r = svg_i["width"]*MhtnDefaults[
                                "Mhtn_Display_SVG_Values"]["GraphValues"]["min_circle_frac"];

    let graph_info = {
        "org": [x_org, y_org], 
        "svg_d3": svg_d3,
        "min_circle_radius": min_crc_r,
        "max_circle_radius": max_crc_r,
        "x_axis_len": x_axis_len,
        "y_axis_len": y_axis_len,
        "max_y_val": y_ticks[y_ticks.length - 1],
        "total_bp_len":  current_scf_info_d["scaffold_length"],
        "scf_to_color": scf_to_color_obj 
    };

    let current_scf_names_and_lengths = [[scf_name, current_scf_info_d["scaffold_length"], 0]]; 
    PopulateManhattanSkyline(current_scf_names_and_lengths, 
            MH_data["scaffolds"], graph_info,
            false, MH_data); 

}



function GetProperTicks(start_val, end_val) {
    /*
    This function is to get properly spread ticks between
    two values, primarily on the y axis.
    start_val: Number
    end_val Number

    Returns:
        ticks_list = [start_val, start_val + subdivs, start_val + 2subdivs,..., end_val]
    */
    subdivs = ConvertValueIntoSubDivs(end_val - start_val);
    tick_values = GetTickValues(start_val, end_val, subdivs);
    return tick_values;
}


function ConvertValueIntoSubDivs(Val) {
    /*
    Important Questions:
    1. Max ticks in axis assuming no more than 3 digits per value?
        Answer: 16
    2. Min ticks in axis?
        Answer: 8

    Meaning: 
        if N = d * 10^n, d > 5 implies division is 5 * 10^(n-2)
        4 < d < 5 implies division is  2.5 * 10^(n-2)
        2 < d < 4 implies division is  2 * 10^(n-2)
        1 < d < 2 implies division is 1 * 10^(n-2)
    */

    let val_info = BaseNotation(Val, 10, 20);
    let dig = val_info[0];
    let power = val_info[1];
    let subdivs = null;

    if (power === 0) {
        subdivs = 1 ;
    } else {
            if (dig >=8) { 
            subdivs =  Math.pow(10,power);
            } else if (dig >= 6) { 
            subdivs = 5 * Math.pow(10, power-1);
            } else {
            subdivs = Math.floor(dig) * Math.pow(10, power-1);
            }
    }
    return subdivs;
}



function GetTickValues(start_val, end_val, subdivs) {

    /*We go from a value and subdivs to actual graph ticks


    Args:
        start_val: (int)
        end_val: (int)
        subdivs: (int)

    Returns:
        ticks_list = [start_val, start_val + subdivs, start_val + 2subdivs,...]

    Specifically, this function starts from start_val and adds subdiv until reaching
        end_val. Note that difference between start_val and end_val does not 
        need t
    */
    // First we get a list of just each tick, not the start and end ticks (no dbl)
    let init_tick_list = [start_val];

    let crnt_val = start_val + subdivs;

    while (crnt_val < end_val){
        init_tick_list.push(crnt_val);
        crnt_val = crnt_val + subdivs;
    }

    init_tick_list.push(end_val);


    return init_tick_list;

}



function BaseNotation(N, base, max_power) {

    /* We get power of base and digit multiplier.
        Eg. if N = 346, base = 10 we return [3.46, 2] because
            3.46 * 10^2 = 346 


    Args:
        N: int, number to find bounds for. MUST BE > 0
        base: int 
        max_power: int (limit so function doesn't run forever with while loop)

    Returns:
        [a, b (power of 10)] where a*10^b = N
        OR [-1, -1] if it failed for some reason

    */

    if (N <= 0) {
        return [-1, -1]
    }
    for (i=0; i < max_power + 1 ;i++){
        if (N >= Math.pow(base,i) && N < Math.pow(base,i+1)) {
            return [ N/Math.pow(base,i), i]
        }
    }

    return [-1, -1]

}

//END CODE IN USE




function BrushCode() {

        //BRUSH CODE: CURRENTLY NOT IN USE

         // Add a clipPath: everything out of this area won't be drawn.
        var clip = svg_obj.append("defs").append("svg:clipPath")
            .attr("id", "clip")
            .append("svg:rect")
            .attr("width", svg_i["width"] )
            .attr("height", svg_i["height"] )
            .attr("x", 0)
            .attr("y", 0);
         

        // Add brushing - making global var
         brush = d3.brushX()                 // Add the brush feature using the d3.brush function
            .extent( [ [0,0], [svg_i["width"], svg_i["height"]] ] ) 
        // initialise the brush area: start at 0,0 and finishes at width,height: it means I select the whole graph area
            .on("end", updateChart) // Each time the brush selection changes, trigger the 'updateChart' function


        // Color scale: give me a specie name, I return a color
        var color = d3.scaleOrdinal()
        .domain(["setosa", "versicolor", "virginica" ])
        .range([ "#440154ff", "#21908dff", "#fde725ff"])

    // Create the scatter variable: where both the circles and the brush take place (global)
    scatter = svg_obj.append('g')
    .attr("clip-path", "url(#clip)")

      // Add the brushing
        scatter
            .append("g")
            .attr("class", "brush")
            .call(brush);

}






function prep_int(inp_i){
    /*
     * Converts an integer input into a comma separated string
     * i.e. 1000 -> 1,000
     * Currently Python
     *
     * Args:
     *     inp_i: Input integer 
     *
     * 
     */

    
    int_l = inp_i.toString().split('');

    op_str = '';
    while (int_l.length > 3) {
        c_char = int_l.shift();
        op_str += c_char;
        if (int_l.length % 3 == 0) {
            op_str += ",";
        }
    }

    op_str += int_l.join('') 

    return op_str
}


/*
 *
 *
 *
function MakeExpandScaffoldButton(svg_obj, tlc_x, width, tlc_y, height, scaffold_name) {
    /*
     *  Q: Are x and y top left corners?
     *
     */
/*
            svg_obj.append('rect')
                .attr('x', tlc_x)
                .attr('y', tlc_y)
                .attr('width', width)
                .attr('height', height)
                .attr('fill', "paleturquoise")
                .attr("cursor", "pointer")
                .data([scaffold_name])
                .on("click", function(d) {
                    MakeSingleScaffoldGraph(svg_obj, d)
                });

            // We create plus symbol - Vertical line
            let c_l = MakeLine(svg_obj, "black", 
                    tlc_x + (width/2), 
                    tlc_y + 2, 
                    tlc_x + (width/2), 
                    tlc_y + height - 2, 
                    .5 )

            
            c_l.attr("cursor", "pointer")
                .data([scaffold_name])
                .on("click", function(d) {
                    MakeSingleScaffoldGraph(svg_obj, d)
                });
                

            // Plus symbol - Horizontal line
            let n_l = MakeLine(svg_obj, "black", 
                    tlc_x + 2, 
                    tlc_y + (height/2), 
                    tlc_x + width - 2, 
                    tlc_y + (height/2), 
                    .5)

            
            n_l.attr("cursor", "pointer")
                .data([scaffold_name])
                .on("click", function(d) {
                    MakeSingleScaffoldGraph(svg_obj, d)
                });
                }
       */     

/*
 *
 *
 */
