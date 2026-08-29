(function () {
    const dropzone = document.getElementById("dropzone");
    const fileInput = document.getElementById("fileInput");
    const fileSelected = document.getElementById("fileSelected");
    const fileName = document.getElementById("fileName");
    const analyzeBtn = document.getElementById("analyzeBtn");
    const loadingStatus = document.getElementById("loadingStatus");
    const errorStatus = document.getElementById("errorStatus");
    const results = document.getElementById("results");

    dropzone.addEventListener("click", () => fileInput.click());

    dropzone.addEventListener("dragover", (event) => {
        event.preventDefault();
        dropzone.classList.add("dragover");
    });

    dropzone.addEventListener("dragleave", () => {
        dropzone.classList.remove("dragover");
    });

    dropzone.addEventListener("drop", (event) => {
        event.preventDefault();
        dropzone.classList.remove("dragover");
        if (event.dataTransfer.files.length) {
            fileInput.files = event.dataTransfer.files;
            updateFileSelection();
        }
    });

    fileInput.addEventListener("change", updateFileSelection);
    analyzeBtn.addEventListener("click", analyzeFile);

    function updateFileSelection() {
        if (fileInput.files.length) {
            fileName.textContent = fileInput.files[0].name;
            fileSelected.hidden = false;
        } else {
            fileSelected.hidden = true;
        }
    }

    function setBusy(busy) {
        analyzeBtn.disabled = busy;
        loadingStatus.hidden = !busy;
        if (busy) {
            errorStatus.hidden = true;
            results.style.display = "none";
        }
    }

    function showError(message) {
        errorStatus.textContent = message;
        errorStatus.hidden = false;
        results.style.display = "none";
    }

    function riskClass(level) {
        const normalized = (level || "").toLowerCase().replace(/\s+/g, "-");
        if (normalized === "high") return "risk-high";
        if (normalized === "moderate") return "risk-moderate";
        if (normalized === "low" || normalized === "potentially-effective") return "risk-low";
        return "risk-moderate";
    }

    function escapeHtml(text) {
        const div = document.createElement("div");
        div.textContent = text;
        return div.innerHTML;
    }

    function analyzeFile() {
        if (!fileInput.files.length) {
            showError("Select a FASTA file before running the analysis.");
            return;
        }

        const formData = new FormData();
        formData.append("file", fileInput.files[0]);
        setBusy(true);

        fetch("/upload", { method: "POST", body: formData })
            .then(async (response) => {
                const data = await response.json().catch(() => ({}));
                setBusy(false);
                if (!response.ok) {
                    showError(data.error || "Request failed (" + response.status + ").");
                    return;
                }
                if (data.error) {
                    showError(data.error);
                    return;
                }
                displayResults(data);
            })
            .catch((err) => {
                setBusy(false);
                showError("Network error: " + err.message);
            });
    }

    function displayResults(data) {
        errorStatus.hidden = true;

        const similarity = data.similarity_score || 0;
        if (similarity < 70) {
            showError(
                "Sequence similarity to SARS-CoV-2 reference is " +
                    similarity +
                    "% (minimum 70%). Upload a SARS-CoV-2 genome FASTA."
            );
            return;
        }

        document.getElementById("validationNote").textContent =
            "SARS-CoV-2 sequence confirmed (" + similarity + "% reference similarity).";

        document.getElementById("seqLength").textContent =
            (data.sequence_info?.length || 0).toLocaleString() + " bp";
        document.getElementById("gcContent").textContent =
            (data.sequence_info?.gc_content || 0) + "%";

        const variant = data.variant?.name
            ? data.variant.name + " (" + (data.variant.confidence || "N/A") + ")"
            : "Not classified";
        document.getElementById("variantName").textContent = variant;

        const mutations = data.mutations?.list || [];
        document.getElementById("mutationList").textContent = mutations.length
            ? mutations.join(", ")
            : "None detected in monitored spike sites";

        const warningEl = document.getElementById("analysisWarning");
        if (data.warning) {
            warningEl.textContent = data.warning;
            warningEl.hidden = false;
        } else {
            warningEl.hidden = true;
        }

        if (data.neutralization_results) {
            const nr = data.neutralization_results;
            document.getElementById("ic50Value").textContent =
                nr.cocktail_prediction + " \u03bcg/mL";

            const riskEl = document.getElementById("escapeRisk");
            riskEl.textContent = nr.summary_risk;
            riskEl.className = "metric-value " + riskClass(nr.summary_risk);

            const tbody = document.getElementById("antibodyRows");
            tbody.innerHTML = "";

            Object.entries(nr.individual_analysis || {}).forEach(([name, info]) => {
                const row = document.createElement("tr");
                row.innerHTML =
                    "<td>" + escapeHtml(name) + "</td>" +
                    "<td>" + escapeHtml(String(info.predicted_ic50_ug_ml ?? "N/A")) + "</td>" +
                    "<td>" + escapeHtml(String(info.fold_change ?? "N/A")) + "</td>" +
                    '<td class="' + riskClass(info.risk_level) + '">' +
                    escapeHtml(info.risk_level || "N/A") +
                    "</td>";
                tbody.appendChild(row);
            });
        }

        results.style.display = "block";
    }
})();
