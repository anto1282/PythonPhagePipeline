

process SAMPLELISTMAKER{
    input:
    val fromVal
    val toVal
    val byVal

    output:
    stdout

    script:
    """
    python3 SampleList.py ${fromVal} ${toVal} ${byVal}
    """
}