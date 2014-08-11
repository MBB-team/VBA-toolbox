---
layout: page
title: Get the toolbox
---

<h2>Download</h2>

<div class="download-option">
        <div class="button">
            <a class="fa fa-archive fa-5x" href="https://github.com/MBB-team/VBA-toolbox/archive/master.zip"></a>
        </div>
        <div class="description">
        <h1>Archive</h1>
        <p>Downlad the last stable release as a .zip file</p>
        </div>
</div>

<div class="download-option">
        <div class="button">
            <a class="fa fa-github-alt fa-5x" href="https://github.com/MBB-team/VBA-toolbox"></a>
        </div>
        <div class="description">
        <h1>Github</h1>
        <p>Clone the repo and stay on the edge</p>
        </div>
</div>


<h2>Install</h2>
<p> Add the VBA-toolbox directory to your Matlab path and save it for automatic loading:
{% highlight matlab %}
addpath(genpath('path_to_the/VBA-toolbox'));
 savepath();
{% endhighlight %}
</p>

<h2>Try it!</h2>
<p> You re ready to go! Launch the <code>demo_Qlearning</code> for a test, and visit the <a href="/wiki">wiki</a> for more informations. </p>



