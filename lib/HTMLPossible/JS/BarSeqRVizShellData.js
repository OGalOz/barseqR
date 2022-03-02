window.BasicShellData = {
    "lyt_vls": {
        "largest_container": {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "body",
                "id": "largest-container-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0.05,
                "t": 0.05,
                "h": 0.9,
                "w": 0.9
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray"
            }
        },
        "full_graph_div": {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "largest-container-div",
                "id": "graph-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": .85,
                "w": 0.7
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray"
            }
        },
        "graph_title_div": {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "graph-div",
                "id": "graph-title-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": .1,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "fontWeight": "bold",
                "textAlign": "center",
                "justifyContent": "center",
                "display": "flex",
                "alignItems": "center",
                "textDecoration": "underline",
                "fontSize": "30px"

            },
            "unq_prp": {
                "innerHTML": "Fitness vs T Score Volcano Plot"
            }
        },
        "graph_svg_container_div": {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "graph-div",
                "id": "graph-svg-container-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": .1,
                "h": .9,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",

            }
        },
        "graph_svg": {
            "tag_type": "SVG",
            "id_i": {
                "parent_id": "graph-svg-container-div",
                "id": "graph-svg",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": 1,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
            }
        },
        "sidebar": {
            "tag_type": "div",
            "id_i": {
                "parent_id": "largest-container-div",
                "id": "graph-sidebar",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0.7,
                "t": 0,
                "h": 1,
                "w": 0.3
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "textAlign": "center"
            }
        },
        "volcano_button_div": {
            "tag_type": "div",
            "id_i": {
                "parent_id": "graph-sidebar",
                "id": "volcano-btn-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": 0.1,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "fontWeight": "bold",
                "textAlign": "center",
                "justifyContent": "center",
                "display": "flex",
                "alignItems": "center",
                "textDecoration": "underline",
                "fontSize": "30px",
                "backgroundColor": "#7fffd4",
                "cursor": "pointer"
            },
            "unq_prp": {
                "innerHTML": "Volcano Plot"
            }
        },
        "compare_conditions_button":{
            "tag_type": "div",
            "id_i": {
                "parent_id": "graph-sidebar",
                "id": "compare-btn-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0.1,
                "h": 0.1,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "fontWeight": "bold",
                "textAlign": "center",
                "justifyContent": "center",
                "display": "flex",
                "alignItems": "center",
                "textDecoration": "underline",
                "fontSize": "30px",
                "cursor": "pointer",
                "backgroundColor": "#00ffff"
            },
            "unq_prp": {
                "innerHTML": "Compare Conditions"
            }
        },
        "chosen_conditions_container": {
            "tag_type": "div",
            "id_i": {
                "parent_id": "graph-sidebar",
                "id": "chosen-full-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0.2,
                "h": 0.2,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
            }
        },
        "chosen_title":{
            "tag_type": "div",
            "id_i": {
                "parent_id": "chosen-full-div",
                "id": "chosen-title-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": 0.2,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "fontWeight": "bold",
                "textAlign": "center",
                "justifyContent": "center",
                "display": "flex",
                "alignItems": "center",
                "textDecoration": "underline",
                "fontSize": "20px"
            },
            "unq_prp": {
                "innerHTML": "Selected Conditions:"
            }
        },
        "chosen_condition_display": {
            "tag_type": "div",
            "id_i": {
                "parent_id": "chosen-full-div",
                "id": "chosen-display-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0.2,
                "h": 0.8,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "overflow": "auto",
                "fontWeight": "normal",
                "fontSize": "15px"
            }
        },
        "clear_selected_btn_div": {
            "tag_type": "div",
            "id_i": {
                "parent_id": "graph-sidebar",
                "id": "clear-selected-btn-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0.4,
                "h": 0.05,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "fontWeight": "bold",
                "cursor": "pointer",
                "textAlign": "center",
                "justifyContent": "center",
                "display": "flex",
                "alignItems": "center",
                "textDecoration": "underline",
                "backgroundColor": "#f0ebfa", 
                "fontSize": "20px"
            },
            "unq_prp": {
                "innerHTML": "Clear Selected Condition"
            }
        },

        "experiments_title_div": {
            "tag_type": "div",
            "id_i": {
                "parent_id": "graph-sidebar",
                "id": "exp-title-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0.45,
                "h": 0.05,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "fontWeight": "bold",
                "textAlign": "center",
                "justifyContent": "center",
                "display": "flex",
                "alignItems": "center",
                "textDecoration": "underline",
                "fontSize": "30px"
            },
            "unq_prp": {
                "innerHTML": "Conditions:"
            }
        },
        "experiments_table_div": {
            "tag_type": "div",
            "id_i": {
                "parent_id": "graph-sidebar",
                "id": "exp-table-div",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0.5,
                "h": 0.5,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "overflow": "auto"
            }
        },

        "experiments_table": {
            "tag_type": "table",
            "id_i": {
                "parent_id": "exp-table-div",
                "id": "exp-table",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": 1,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "overflow": "auto"
            }
        },
        "bottom_info_bar": {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "largest-container-div",
                "id": "bottom-bar",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": .85,
                "h": .15,
                "w": 0.7
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray"
            }
        },
        "selected_point_info_cnt" : {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "bottom-bar",
                "id": "selected-point-full",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": 1,
                "w": 0.3
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray"
            }
        },
        "selected_point_info_title" : {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "selected-point-full",
                "id": "selected-point-title",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0,
                "h": .3,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "fontWeight": "bold",
                "textAlign": "center",
                "justifyContent": "center",
                "display": "flex",
                "alignItems": "center",
                "textDecoration": "underline",
                "fontSize": "18px"
            },
            "unq_prp": {
                "innerHTML": "Selected point info:"
            }
        },
        "selected_point_info_display" : {
            "tag_type": "DIV",
            "id_i": {
                "parent_id": "selected-point-full",
                "id": "selected-point-display",
                "class": "base-shell"
            },
            "size_loc_i": {
                "values_type": "fractions",
                "l": 0,
                "t": 0.3,
                "h": 0.7,
                "w": 1
            },
            "style_i": {
                "position": "absolute",
                "border": "2px solid gray",
                "overflow": "auto"
            }
        }


        

    }
}


