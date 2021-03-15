/* Len Taing 2020 (TGBTG) */

/* Function to toggle sections in the main panel */
//var last_section = 'meta';
function toggler(divId) {
    if (last_section != null) {
        $("#" + last_section).toggle();
    }
    $("#" + divId).toggle();
    last_section = divId;
}

//HANDLE table cells image clicks
$('.wes-image-modal').on('click',function(){
    $('#wesImageModal').modal({show:true});
    $('#wesImageModal_img').attr("src", this.src);
});

//Create datatables
$('.wes_datatable').each(function(index) {
    //Create the DataTable table
    var originalTbl = $(this).clone(); //copy used in downloadBtn hdlr below
    var dt = $(this).DataTable({lengthMenu: [[10, 50, 100, -1], [10, 50, 100, "All"]]});
    var tblTitle = $(this).attr('id');

    //Insert a Download list btn JUST after the 'Show entries' menu
    var entryPt = $(this).parent().find('.dataTables_length');
    var downloadBtn = $('<button class=\"btn btn-primary btn-sm ml-3\">Download CSV</button>');
    entryPt.append(downloadBtn);

    //Copy table contents and generate new DataTable in hidden modal
    //invoke the download event
    downloadBtn.on('click', function() {
	//Clone table
	var newTbl = originalTbl.clone();
	//Add it to hiddent modal
	mb = $('#wesSubModal').find('.modal-body');
	mb.empty();
	mb.append(newTbl);
	//Create new DT
	var tmp = newTbl.DataTable({dom: 'Bfrtip', buttons: [{extend:'csvHtml5', title: tblTitle}]});
	//invoke download btn
	tmp.buttons().trigger('click');
    });
});

