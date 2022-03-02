 /* In this file we create a Scatter Plot
     CURRENTLY: drawDashedZeroLines
     */


async function refreshScatterPlotAxes(svg_id, scatter_data, axes_info,
                                     dashed_lines_bool=true) {

    /* Args:
      svg_id: (str) DOM id of svg object
      scatter_data: (Object)
         min_y: Number
         max_y: Number
         min_x: Number 
         max_x: Number
         [point_list]: list<<fit (Num), t (Num)>, ... >
      axes_info: (Object) 
        Like svg_axes_info as defined in 
            makeSVGAxesAndLabels from d3SVGUtil.js
            but with some parts left undefined,
            like x axis left corner
    */

    clearSVG(svg_id)
    new_ret_d = await createScatterPlotAxes(svg_id, scatter_data, axes_info,
                                     dashed_lines_bool)

    window.ret_d = new_ret_d
    return new_ret_d
}


async function createScatterPlotAxes(svg_id, scatter_data, axes_info,
                                     dashed_lines_bool=true) {
    /* Args:
      svg_id: (str) DOM id of svg object
      scatter_data: (Object)
         min_y: Number
         max_y: Number
         min_x: Number 
         max_x: Number
         [point_list]: list<<fit (Num), t (Num)>, ... >
      axes_info: (Object) 
        Like svg_axes_info as defined in 
            makeSVGAxesAndLabels from d3SVGUtil.js
            but with some parts left undefined,
            like x axis left corner
    */

    let x_axis_start_val = Math.ceil(scatter_data["min_x"] - 1);
    let x_axis_end_val = Math.floor(scatter_data["max_x"] + 1);
    let y_axis_start_val = Math.floor(scatter_data["min_y"] - 1);
    let y_axis_end_val = Math.ceil(scatter_data["max_y"] + 1);

    // GetProperTicks is function from d3SVGUtil
    let x_ticks_list = GetProperTicks(x_axis_start_val, x_axis_end_val)
    let y_ticks_list = GetProperTicks(y_axis_start_val, y_axis_end_val)
    // ^ Note x_ticks_list and y_ticks_list are list<Num>

    // We add the values 0 to both the x_ticks_list and the y_ticks_list
    let x_info_a = addZeroToTicksList(x_ticks_list)
    let y_info_a = addZeroToTicksList(y_ticks_list)
    // We get ticks lists and whether or not they cross zero
    x_ticks_list = x_info_a[0]
    let x_crosses_zero = x_info_a[1]
    y_ticks_list = y_info_a[0]
    let y_crosses_zero = y_info_a[1]
    
    // Just the main axes and labels for the axes, no ticks
    let return_object = makeSVGAxesAndLabels(svg_id, axes_info)

    /*
     return_object contains:
        x_axis_len: Num
        y_axis_len: Num
        graph_blc: (coordinates [x,y])
        graph_trc: (coordinates [x,y])
        y_axis_start: (coordinates [x,y])
        x_axis_start: (coordinates [x,y])
    */

    // We draw the ticks and the labels for X, then Y:

    // First we need to generate the tick_info_list:
    
    x_tick_info_list = makeStandardTickInfoListFromTicks(x_ticks_list, 10)
    makeAxisTicksAndTickLabels(svg_id, "x",
                               return_object['x_axis_start'], 
                               [return_object['x_axis_start'][0] + return_object['x_axis_len'],
                               return_object['x_axis_start'][1]],
                               x_tick_info_list) 

    y_tick_info_list = makeStandardTickInfoListFromTicks(y_ticks_list, 10)
    makeAxisTicksAndTickLabels(svg_id, "y",
                               return_object['y_axis_start'], 
                               [return_object['y_axis_start'][0],
                               return_object['y_axis_start'][1]  - return_object['y_axis_len']],
                               y_tick_info_list) 


    if (dashed_lines_bool) {
        // We draw dashed lines where the 0's of each axis are
        d3svg = getSVG(svg_id);
        drawDashedZeroLines(d3svg, x_ticks_list, y_ticks_list,
                            return_object['x_axis_len'], return_object['x_axis_start'],
                            return_object['y_axis_len'], return_object['y_axis_start'],
                            x_crosses_zero, y_crosses_zero)
    }

    return {
        "x_ticks_list" : x_ticks_list,
        "y_ticks_list" : y_ticks_list,
        "x_axis_len" : return_object['x_axis_len'],
        "y_axis_len" : return_object['y_axis_len'],
        "x_axis_start" : return_object['x_axis_start'],
        "y_axis_start" : return_object['y_axis_start'],
        "graph_blc": return_object['graph_blc'],
        "graph_trc": return_object['graph_trc']
    }

}



async function populateSVGGraphWithPoints(svg_id,
                                          point_list,
                                          x_ticks_list,
                                          y_ticks_list,
                                          x_axis_start,
                                          y_axis_start,
                                          x_axis_len,
                                          y_axis_len,
                                          async_points_val=100000) {

    printLoadingSign('bottom-bar', 0)
    /*
    let ld = setInterval( function() {
        printLoadingSign('bottom-bar', 1)
    }, 500)
    */

    // async part of function:
    setTimeout(function() {
       let finished = delayedPopulateScatterPlot(d3svg,
                        point_list,
                        x_ticks_list,
                        y_ticks_list,
                        ret_d['x_axis_start'],
                        ret_d['y_axis_start'],
                        ret_d['x_axis_len'],
                        ret_d['y_axis_len'],
                        300000)
       finished.then((value) => {
           console.log("Finished plotting") 
           // We show the user we are done loading points. Remove sign
           // (from d3SVGUtil)
           printLoadingSign('bottom-bar', 2)
           //window.still_loading = false
           //clearInterval(ld)
       })

    }, 10)

}


async function delayedPopulateScatterPlot(d3svg, 
                                      point_list,
                                      x_ticks_list,
                                      y_ticks_list,
                                      x_axis_start,
                                      y_axis_start,
                                      x_axis_len,
                                      y_axis_len,
                                      split_num=50000) {

    // Returns a promise if it has few enough points to update 

    pl_len = point_list.length; 
    if (pl_len > 1000000000) {
        console.log(pl_len)
        throw "point_list length is too long, should be less than 10^9"
    } else if (pl_len > split_num) {
        let remainder = pl_len % split_num
        let num_divs = (pl_len - remainder)/split_num
        console.log(num_divs)
        window.current_div = -1
        window.new_current = 0
        var delayedPopulateInterval = setInterval(function() {
                if (window.new_current > window.current_div) {
                    window.current_div = window.new_current;
                    let new_promise = setTimeout(function() {
                        indexedPopulate(current_div,
                                        num_divs,
                                        d3svg,
                                        point_list,
                                        split_num,
                                        x_ticks_list,
                                        y_ticks_list,
                                        x_axis_start,
                                        y_axis_start,
                                        x_axis_len,
                                        y_axis_len)
                    new_promise.then((value) => {
                        window.new_current += 1
                    } )
                    }, 1000)
                } else if (window.new_current == num_divs) {

                    // Goes from last remainder to the end
                    console.log("Plotting " + num_divs.toString() + "/" + num_divs.toString() 
                            + " fractions of the total points")
                    let final_populated = populateSVGWithScatterPoints(d3svg,
                                     point_list.slice(
                                         (num_divs)*split_num
                                     ),
                                     x_ticks_list,
                                     y_ticks_list,
                                     x_axis_start,
                                     y_axis_start,
                                     x_axis_len,
                                     y_axis_len)
                    final_populated.then((value) => {

                        clearInterval(delayedPopulateInterval)
                        printLoadingSign('bottom-bar', 2)
                        return new Promise((resolve, reject) => {
                            resolve(true)
                        })
                    })
                }
            
           printLoadingSign('bottom-bar', 1)
        }, 500)
    
        /*
     window.final_populated.then((value) => {
                return new Promise((resolve, reject) => {
                        resolve(true)
                    })

    })
    */

    } else {
        let populated = populateSVGWithScatterPoints(d3svg,
                             point_list,
                             x_ticks_list,
                             y_ticks_list,
                             x_axis_start,
                             y_axis_start,
                             x_axis_len,
                             y_axis_len
                             )

        populated.then((value) => {
            return new Promise((resolve, reject) => {
                    resolve(true)
                })

        })
    }

}

function indexedPopulate(i,
                max_num,
                d3svg,
                point_list,
                split_num,
                x_ticks_list,
                y_ticks_list,
                x_axis_start,
                y_axis_start,
                x_axis_len,
                y_axis_len) {
    /*
     *
     * Args:
     *  i: the current index of work to be done
     *  max_num: The maximum index of work to be done
     */
    
    if ( i < max_num) {
        
         let populated = populateSVGWithScatterPoints(d3svg,
                                     point_list.slice(
                                         i*split_num,
                                         (i+1)*split_num
                                     ),
                                     x_ticks_list,
                                     y_ticks_list,
                                     x_axis_start,
                                     y_axis_start,
                                     x_axis_len,
                                     y_axis_len
                                 )
        /*
        populated.then((value) => {
            console.log("fulfilled")
        })
        */
        return new Promise((resolve, reject) => {
            resolve(true)
        })
    } else if (i == max_num){
        return new Promise((resolve, reject) => {
            resolve(true)
        })
    }else {
        console.log("i surpassed max_num")

        clearInterval(delayedPopulateInterval)
    }

}




async function populateSVGWithScatterPoints(svg_id,
                                      point_list,
                                      x_ticks,
                                      y_ticks,
                                      x_axis_start,
                                      y_axis_start,
                                      x_axis_len,
                                      y_axis_len,
                                      point_contains_data = false, 
                                      point_click_function = null
                                      ) {
    /*
     *
     * point_list: list<[x (Num), y (Num)]> for all the points
     *      based on actual numerical data, not coordinates
     *  x_ticks: list<Num> ticks in x axis
     *  y_ticks: list<Num> ticks in y axis
     *  x_axis_start: coordinates [Num, Num]
     *  y_axis_start: coordinates [Num, Num]
     *  x_axis_len: Num
     *  y_axis_len: Num
     *
     */

    let d3svg = getSVG(svg_id)

    let total_x_range = x_ticks[x_ticks.length -1] - x_ticks[0];
    let total_y_range = y_ticks[y_ticks.length -1] - y_ticks[0];

    //quadrant coloration
    var my_qc = {
        'q1': 'green',
        'q3': 'purple',
        'q2': 'red',
        'q4': 'blue'
    }
    
    // Note: addManyPointsToPlot is function from d3SVGUtil
    // points_added is a Promise
    let points_added = addManyPointsToPlot(d3svg,
                            point_list,
                            total_x_range,
                            x_axis_len,
                            x_axis_start,
                            x_ticks[0],
                            total_y_range,
                            y_axis_len,
                            y_axis_start,
                            y_ticks[0],
                            3,
                            async=true,
                            quadrant_coloration = my_qc,
                            point_contains_data= point_contains_data,
                            point_click_function= point_click_function
                            )

    points_added.then((value) => {
        console.log("Added points to plot A")
        return new Promise((resolve, reject) => {
            resolve(true)
        })
    })

}




function nullfunc(x) {
    //null func
    
}


function drawDashedZeroLines(d3svg, x_ticks_list, y_ticks_list, 
                            x_axis_len, x_axis_start,
                            y_axis_len, y_axis_start,
                            draw_x_zero_dashed=true,
                            draw_y_zero_dashed=true) {
    /*
     * We find where 0's occur within the ticks list
     * and their ratio to the total length of the axes
     * and draw a dashed line.
     * To find where the 0's occur, we take the full
     * numerical distance between the tick's lowest value
     * and the ticks highest value and find range of 0 from
     * lowest value and create line perpendicularly from
     * that point.
     *  
     *
        x_axis_len: Num
        y_axis_len: Num
        y_axis_start: (coordinates [x,y])
        x_axis_start: (coordinates [x,y])
     * We take the function makeDashedLine(d3svg, start_coordinates, end_coordinates,
                        dash_length, break_length,
                        color, stroke_width)
     */


    if (draw_x_zero_dashed) {

        let total_y_range = y_ticks_list[y_ticks_list.length -1] - y_ticks_list[0];
        let x_zero_frac = ((-1)*y_ticks_list[0])/total_y_range;
        let x_zero_dashed_start = [y_axis_start[0], 
                                y_axis_start[1] - y_axis_len*x_zero_frac]
        let x_zero_dashed_end = [y_axis_start[0] + x_axis_len, 
                                y_axis_start[1] - y_axis_len*x_zero_frac]


        makeDashedLine(d3svg, x_zero_dashed_start, x_zero_dashed_end,
                        3,3,
                        "green", 1)
    }

    if (draw_y_zero_dashed) {
        
        let total_x_range = x_ticks_list[x_ticks_list.length -1] - x_ticks_list[0];
        let y_zero_frac = ((-1)*x_ticks_list[0])/total_x_range

        let y_zero_dashed_start = [x_axis_start[0] + x_axis_len*y_zero_frac, 
                                x_axis_start[1]]
        let y_zero_dashed_end = [x_axis_start[0] + x_axis_len*y_zero_frac, 
                                x_axis_start[1] - y_axis_len]

        makeDashedLine(d3svg, y_zero_dashed_start, y_zero_dashed_end,
                        3,3,
                        "blue", 1)
    }




}

function addZeroToTicksList(ticks_list) {
    /* ticks_list is list<Num> to which we add the number 0
        if it doesn't exist
        If the ticks_list doesn't cross zero
        we return crosses_zero_bool = False
        
        Returns:
            <new_tick_list, crosses_zero_bool>
    */

    let transitional_index = -1;
    let less_than_zero = true;
    let crosses_zero_bool = true;

    for (let i=0; i<ticks_list.length; i++) {
        let c_t = ticks_list[i];
        if (c_t >= 0) {
            if (c_t == 0) {
                return [ticks_list, crosses_zero_bool];
            } else {
                    transitional_index = i;
                    less_than_zero = false;
                    break;
            }
        }
    }
    
    if (less_than_zero == true) {
        console.log("ticks_list never crosses zero")
        crosses_zero_bool = false;
        return [ticks_list, crosses_zero_bool]
    } else {
        let new_ticks_list = ticks_list.slice(0,transitional_index).concat([0]) 
                                    .concat(ticks_list.slice(transitional_index));
        return [new_ticks_list, crosses_zero_bool]
    }
    
}


/*
async function populateSVGWithScatterPlot(svg_id, scatter_data, axes_info) {
    /* Args:
      svg_id: (str) DOM id of svg object
      scatter_data: (Object)
         min_y: Number
         max_y: Number
         min_x: Number 
         max_x: Number
         point_list: list<<fit (Num), t (Num)>, ... >
      axes_info: (Object) 
        Like svg_axes_info as defined in 
            makeSVGAxesAndLabels from d3SVGUtil.js
            but with some parts left undefined,
            like x axis left corner
    */
/*
    
    let x_axis_start_val = Math.ceil(scatter_data["min_x"] - 1);
    let x_axis_end_val = Math.floor(scatter_data["max_x"] + 1);
    let y_axis_start_val = Math.floor(scatter_data["min_y"] - 1);
    let y_axis_end_val = Math.ceil(scatter_data["max_y"] + 1);

    // GetProperTicks is function from d3SVGUtil
    let x_ticks_list = GetProperTicks(x_axis_start_val, x_axis_end_val)
    let y_ticks_list = GetProperTicks(y_axis_start_val, y_axis_end_val)
    // ^ Note x_ticks_list and y_ticks_list are list<Num>

    // We add the values 0 to both the x_ticks_list and the y_ticks_list
    let x_info_a = addZeroToTicksList(x_ticks_list)
    let y_info_a = addZeroToTicksList(y_ticks_list)
    // We get ticks lists and whether or not they cross zero
    x_ticks_list = x_info_a[0]
    // below is a bool
    let x_crosses_zero = x_info_a[1]
    y_ticks_list = y_info_a[0]
    // below is a bool
    let y_crosses_zero = y_info_a[1]
    
    // Just the main axes and labels for the axes, no ticks
    ret_d = makeSVGAxesAndLabels(svg_id, axes_info)


    /*
     ret_d contains:
        x_axis_len: Num
        y_axis_len: Num
        graph_blc: (coordinates [x,y])
        graph_trc: (coordinates [x,y])
        y_axis_start: (coordinates [x,y])
        x_axis_start: (coordinates [x,y])
    */
/*
    // We draw dashed lines where the 0's of each axis are
    d3svg = getSVG(svg_id);
    drawDashedZeroLines(d3svg, x_ticks_list, y_ticks_list,
                        ret_d['x_axis_len'], ret_d['x_axis_start'],
                        ret_d['y_axis_len'], ret_d['y_axis_start'],
                        x_crosses_zero, y_crosses_zero)


    // We draw the ticks and the labels for X, then Y:

    // First we need to generate the tick_info_list:
    
    x_tick_info_list = makeStandardTickInfoListFromTicks(x_ticks_list, 10)
    makeAxisTicksAndTickLabels(svg_id, "x",
                               ret_d['x_axis_start'], 
                               [ret_d['x_axis_start'][0] + ret_d['x_axis_len'],
                               ret_d['x_axis_start'][1]],
                               x_tick_info_list) 

    y_tick_info_list = makeStandardTickInfoListFromTicks(y_ticks_list, 10)
    makeAxisTicksAndTickLabels(svg_id, "y",
                               ret_d['y_axis_start'], 
                               [ret_d['y_axis_start'][0],
                               ret_d['y_axis_start'][1]  - ret_d['y_axis_len']],
                               y_tick_info_list) 



}
*/
