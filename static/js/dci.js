
document.getElementById("file_num").addEventListener("input", function() {
    var num = parseInt(this.value);
    if (num < 2) {
        alert("File num should greater than 1.");
        document.getElementById("file_num").value = 2;
        num = 2;
    }
    var container = document.getElementById("file_inputs_container");
    container.innerHTML = ""; // 清空容器

    var t = document.createElement("table");
    t.id = "file_inputs_table";

    // 创建表头
    var thead = document.createElement("thead");
    var tr = document.createElement("tr");
    var th1 = document.createElement("th");
    th1.textContent = "file";
    var th2 = document.createElement("th");
    th2.textContent = "sample_name";
    var th3 = document.createElement("th");
    th3.textContent = "group";
    tr.appendChild(th1);
    tr.appendChild(th2);
    tr.appendChild(th3);
    thead.appendChild(tr);
    t.appendChild(thead);

    // 创建表格内容
    var tbody = document.createElement("tbody");
    for (var i = 0; i < num; i++) {
        var tr = document.createElement("tr");
        var td1 = document.createElement("td");
        var td2 = document.createElement("td");
        var td3 = document.createElement("td");
        var input1 = document.createElement("input");
        input1.id = "file_input_" + i;
        input1.type = "text";
        var input2 = document.createElement("input");
        input2.id = "sample_name_input_" + i;
        input2.type = "text";
        var input3 = document.createElement("input");
        input3.id = "group_input_" + i;
        input3.type = "text";
        td1.appendChild(input1);
        td2.appendChild(input2);
        td3.appendChild(input3);
        tr.appendChild(td1);
        tr.appendChild(td2);
        tr.appendChild(td3);
        tbody.appendChild(tr);
    }
    t.appendChild(tbody);

    container.appendChild(t);
});

function run_dci() {
    console.log('[FUNC] run_dci')// 创建一个XMLHttpRequest对象

    // a loading gif in /static/images/donut.gif
    result_div = document.getElementById("result");
    window.scrollTo({"behavior": "smooth", "top": result_div.offsetTop});

    var num = document.getElementById("file_num").value;
    var files = [];
    var samples = [];
    var groups = [];
    for (var i = 0; i < num; i++) {
        var file = document.getElementById("file_input_" + i).value;
        var sample = document.getElementById("sample_name_input_" + i).value;
        var group = document.getElementById("group_input_" + i).value;
        files.push(file)
        samples.push(sample)
        groups.push(group)
    }
    if (files.length < 2) {
        alert("Please input files.");
        return;
    }
    else if (files.indexOf("") != -1) {
        alert("Please fill full files.")
        return;
    }
    if (samples.length < 2) {
        alert("Please input samples.");
        return;
    }
    else if (samples.indexOf("") != -1) {
        alert("Please fill full samples.")
        return;
    }
    if (groups.length < 2) {
        alert("Please input groups.");
        return;
    }
    else if (groups.indexOf("") != -1) {
        alert("Please fill full groups.")
        return;
    }
    var blacklist_path = document.getElementById("blacklist_path").value;
    var gap = document.getElementById("gap").value;
    if (gap == "") {
        gap = 500;
    }
    var remove_loops_in_blacklist = document.getElementById("remove_loops_in_blacklist").checked;
    var remove_self_ligated_loops = document.getElementById("remove_self_ligated_loops").checked;
    var fdr_threshold = document.getElementById("fdr_threshold").value;
    if (fdr_threshold == "") {
        fdr_threshold = 0.01;
    }
    var loop_len_threshold = document.getElementById("loop_len_threshold").value;
    if (loop_len_threshold == "") {
        loop_len_threshold = 5000;
    }
    var intra_only = document.getElementById("intra_only").checked;
    var chr_filter = document.getElementById("chr_filter").value;
    if (chr_filter == "") {
        chr_filter = "chrM,chrX";
    }
    var pet_threshold = document.getElementById("pet_threshold").value;
    if (pet_threshold == "") {
        pet_threshold = 5;
    }

    if (document.getElementById("norm_checkbox").checked) {
        var norm_method = document.getElementById("norm_method").value;
        var fig_demand = document.getElementById("fig_demand").checked;
    }
    else {
        var norm_method = false;
        var fig_demand = false;
    }
    var reproducibility_checkbox = document.getElementById("reproducibility_checkbox").checked;
    var idr_checkbox = document.getElementById("idr_checkbox").checked;
    if (document.getElementById("dci_checkbox").checked) {
        var dci_method = document.getElementById("dci_method").value;
    }
    else {
        var dci_method = false;
    }

    var output_path = document.getElementById("output_path").value;

    var data = JSON.stringify({ files: files,
        samples: samples, groups: groups, blacklist_path: blacklist_path,
        gap: gap, remove_loops_in_blacklist: remove_loops_in_blacklist, remove_self_ligated_loops: remove_self_ligated_loops,
        fdr_threshold: fdr_threshold, loop_len_threshold: loop_len_threshold, intra_only: intra_only,
        chr_filter: chr_filter, pet_threshold: pet_threshold,
        fig_demand: fig_demand,
        norm_method: norm_method, reproducibility_checkbox: reproducibility_checkbox, idr_checkbox: idr_checkbox, dci_method: dci_method,
        output_path: output_path });

    result_div.innerHTML = "<img src='/static/images/donut.gif'>";
    var xhr = new XMLHttpRequest();

    // 设置请求方法和URL
    xhr.open('POST', '/run_dci', true);
    xhr.setRequestHeader('Content-Type', 'application/json;charset=UTF-8');

    // 设置请求完成时的回调函数
    xhr.onload = function() {
        if (xhr.status === 200) {
            // 请求成功，处理响应数据
            var response = JSON.parse(xhr.responseText);
            if (response["success"] == true) {
                // result_div.innerHTML = "DCI analysis completed!";
                result_div.innerHTML = "Analysis completed!<br>";
                if (document.getElementById("reproducibility_checkbox").checked) {
                    var reproducibility_result_table = "<table id='reproducibility_result'><tr><th>replicate 1</th><th>replicate 2</th><th>reproducibility</th></tr>";
                    for (var i = 0; i < response["aredci_result"].length; i++) {
                        reproducibility_result_table += "<tr><td>" + response["aredci_result"][i][0] + "</td><td>" + response["aredci_result"][i][1] + "</td><td>" + response["aredci_result"][i][2] + "</td></tr>";
                    }
                    reproducibility_result_table += "</table>";
                    result_div.innerHTML += "<h2>Reproducibility result</h2>" + reproducibility_result_table;
                    result_div.innerHTML += "<hr/><h2>DCI result</h2>";
                }
                else {
                    result_div.innerHTML += "<div>Reproducibility assessment was not executed.</div>";
                }
                if (document.getElementById("dci_checkbox").checked) {
                    dci_result_desc = "<ul><li>Loops filtered by Quality Control: "+output_path+"/loops_filtered_in_QC</li>";
                    dci_result_desc = "<li>Intermediate loops after QC: "+output_path+"/[cellline]_[rep]_replace_anchors_test.MICC</li>";
                    dci_result_desc += "<li>Figures generated in Normalization: "+output_path+"/figures</li>";
                    dci_result_desc += "<li>Interactions after Normalization: "+output_path+"/[cellline]_[rep]_dci_data.narrowPeak</li>";
                    dci_result_desc += "<li>DCI identification result: "+output_path+"/my_dci_result.csv</li></ul><br/>";
                    result_div.innerHTML += dci_result_desc;
                }
                else {
                    result_div.innerHTML += "<div>DCI identification was not executed.</div>";
                }
            }
            else {
                result_div.innerHTML = "DCI analysis failed!";
            }
            window.scrollTo({"behavior": "smooth", "top": result_div.offsetTop});
        } else {
            // 请求失败，处理错误情况
            result_div.innerHTML = "DCI analysis failed!";
            console.error('Request failed. Status: ' + xhr.status);
        }  
    };  

    // 设置请求过程中出现错误时的回调函数
    xhr.onerror = function() {
        result_div.innerHTML = "DCI analysis failed!";
        console.error('An error occurred during the request.');
    };

    // 发送请求
    xhr.send(data);
}

function example() {
    var file_num_element = document.getElementById("file_num");
    file_num_element.value = 4;
    // 创建Event对象  
    var event = new Event('input', {  
        'view': window,  
        'bubbles': true,  
        'cancelable': true  
    });
    // 触发事件  
    file_num_element.dispatchEvent(event);
    // document.getElementById("file_input_0").value = "./data/GM12878_ChIA-PET2_result/replicate1/GM12878.interactions.MICC";
    // document.getElementById("file_input_1").value = "./data/GM12878_ChIA-PET2_result/replicate2/GM12878.interactions.MICC";
    document.getElementById("file_input_0").value = "./data/ChIA-PET2_result/K562_ChIA-PET2_result/replicate1/K562.interactions.MICC";
    document.getElementById("file_input_1").value = "./data/ChIA-PET2_result/K562_ChIA-PET2_result/replicate2/K562.interactions.MICC";
    document.getElementById("file_input_2").value = "./data/ChIA-PET2_result/MCF-7_ChIA-PET2_result/replicate1/MCF-7.interactions.MICC";
    document.getElementById("file_input_3").value = "./data/ChIA-PET2_result/MCF-7_ChIA-PET2_result/replicate2/MCF-7.interactions.MICC";
    // document.getElementById("sample_name_input_0").value = "GM12878_rep1";
    // document.getElementById("sample_name_input_1").value = "GM12878_rep2";
    document.getElementById("sample_name_input_0").value = "K562_rep1";
    document.getElementById("sample_name_input_1").value = "K562_rep2";
    document.getElementById("sample_name_input_2").value = "MCF7_rep1";
    document.getElementById("sample_name_input_3").value = "MCF7_rep2";
    // document.getElementById("group_input_0").value = "GM12878";
    // document.getElementById("group_input_1").value = "GM12878";
    document.getElementById("group_input_0").value = "K562";
    document.getElementById("group_input_1").value = "K562";
    document.getElementById("group_input_2").value = "MCF7";
    document.getElementById("group_input_3").value = "MCF7";
    document.getElementById("blacklist_path").value = "./data/blacklist/hg38-blacklist.v2.bed";
    document.getElementById("gap").value = 500;
    document.getElementById("remove_loops_in_blacklist").checked = true;
    remove_self_ligated_loops = document.getElementById("remove_self_ligated_loops").checked = true;
    fdr_threshold = document.getElementById("fdr_threshold").value = 0.01;
    loop_len_threshold = document.getElementById("loop_len_threshold").value = 5000;
    intra_only = document.getElementById("intra_only").checked = true;
    chr_filter = document.getElementById("chr_filter").value = "chrM,chrX";
    document.getElementById("pet_threshold").value = 5;
    document.getElementById("output_path").value = "./data/my_test_data/";
}


{/* <input type="file" id="fileInput" multiple style="display:none">
    <button onclick="selectFiles()">Select Files</button>
    <ul id="selectedFiles"></ul>
    </div>
    <script>
    function selectFiles() {
        const fileInput = document.getElementById('fileInput');
        fileInput.click();
        fileInput.addEventListener('change', handleFileSelection);
    }

    function handleFileSelection(event) {
        const selectedFiles = event.target.files;
        const selectedFilesList = document.getElementById('selectedFiles');
        selectedFilesList.innerHTML = '';
        for (let i = 0; i < selectedFiles.length; i++) {
            const listItem = document.createElement('li');
            listItem.textContent = selectedFiles[i].path || selectedFiles[i].webkitRelativePath; // 获取文件的完整路径
            console.log(selectedFiles[i].path);
            console.log(URL.createObjectURL(selectedFiles[i]));
            console.log(selectedFiles[i].value);
            selectedFilesList.appendChild(listItem);
        }
    }
</script> */}