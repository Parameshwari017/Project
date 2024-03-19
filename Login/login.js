

function validate()
{
    let name1=document.getElementById("name").value

if(name1===""){
    alert("Enter Your UserName")
    document.form1.name.focus()
    return false;
}

if(!isNaN(name1)){
    alert("Enter The Characters / Alphabet Only")
    document.form1.name.value=""
    document.form1.name.focus()
    return false;
}
}


function myFunction() {
    var x = document.getElementById("pwd");
    if (x.type === "password") {
      x.type = "text";
    } else {
      x.type = "password";
    }
  }
